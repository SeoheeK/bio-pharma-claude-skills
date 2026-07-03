---
name: pk-allometry-ensemble
version: 1.0
description: >
  동물→인간 약동학 파라미터 스케일링(Allometry) 앙상블 스킬.
  ChEMBL API에서 종간 PK 데이터를 수집하고, 10가지 Allometry 방법을 
  앙상블하여 인간 CL/Vd를 예측함.
  UniPlot(종간 Log-Log 시각화) + Rule of Exponents + ML 앙상블 포함.

  다음 요청에 반드시 사용:
  - "동물 PK → 인간 예측", "allometry", "종간 스케일링"
  - "CL BW^0.75", "rule of exponents", "MLP 보정", "brain weight 보정"
  - "전임상 → FIH 스케일링", "IVIVE", "mPBPK 스케일링"
  - "ChEMBL PK 데이터 가져와줘", "생체활성 데이터 수집"
  - "UniPlot 그려줘", "종간 PK 비교 그래프"
  - "rat dog monkey → human CL 예측"

  출력물:
  1. ChEMBL 종간 PK 데이터 (parquet/CSV)
  2. 10종 Allometry 방법별 예측값 + 앙상블 결과
  3. UniPlot (Log-Log 분산형 + 앙상블 범위 시각화)
  4. Snakemake 워크플로우 자동 실행 (workflow-01 + workflow-02)

allometry_methods:
  1_simple:           "CL = a × BW^b (단순 멱함수 회귀)"
  2_fixed_exponent:   "b = 0.75 고정 (이론적 allometry)"
  3_rule_of_exponents: "ROE: b<0.7=단순, 0.7-1.0=MLP, >1.0=BrW"
  4_mlp_correction:   "최대 수명(MLP) 보정 — Mahmood 1996"
  5_brain_weight:     "뇌중량(BrW) 보정 — 고 대사 약물"
  6_two_species:      "rat+dog 2종 외삽 — 최소 2종 필요"
  7_fu_corrected:     "단백결합률(fu) 보정 (FCIM)"
  8_ivive:            "In vitro CLint → Well-stirred 간 모델"
  9_mPBPK:            "간 혈류량 기반 Minimal PBPK"
  10_ml_ensemble:     "RandomForest + GradientBoosting LOO 검증"

key_libraries:
  - chembl_webresource_client  # ChEMBL 공식 Python 클라이언트
  - pharmpy                     # NONMEM 파싱 + 분석
  - scikit-learn                # ML 앙상블
  - plotly                      # Interactive UniPlot
  - snakemake                   # 워크플로우 자동화
---

# Allometry Ensemble Skill

## ChEMBL API 연동 가이드

```python
# 설치
# pip install chembl-webresource-client chembl-downloader

from chembl_webresource_client.new_client import new_client

# 1. 약물 ChEMBL ID 조회
molecule = new_client.molecule
results = molecule.filter(pref_name__iexact="busulfan").only(
    ["molecule_chembl_id", "pref_name", "molecule_properties"]
)
chembl_id = results[0]["molecule_chembl_id"]  # "CHEMBL709"

# 2. 종간 PK 파라미터 조회 (CL, Vd)
activity = new_client.activity
pk_data = activity.filter(
    molecule_chembl_id=chembl_id,
    standard_type__in=["CL", "Vd", "HALF-LIFE"],
    target_organism__in=[
        "Homo sapiens",
        "Rattus norvegicus",
        "Mus musculus",
        "Canis lupus familiaris",
        "Macaca fascicularis"
    ]
).only([
    "standard_type", "standard_value", "standard_units",
    "target_organism", "assay_description"
])

# 3. 분자 특성 조회 (MW, LogP)
props = results[0].get("molecule_properties", {})
mw     = props.get("full_mwt")
logp   = props.get("alogp")
hbd    = props.get("hbd")        # H-bond donors
hba    = props.get("hba")        # H-bond acceptors
tpsa   = props.get("psa")        # Polar surface area
ro5    = props.get("num_ro5_violations")  # Lipinski RO5 위반

# 4. 생체활성 데이터 조회 (IC50, EC50)
bio = activity.filter(
    molecule_chembl_id=chembl_id,
    standard_type__in=["IC50", "EC50", "Ki"],
    pchembl_value__gte=5.0,       # pIC50 ≥ 5 (IC50 ≤ 10 µM)
    assay_type="B",               # Binding assay
    target_organism="Homo sapiens"
)
```

## Allometry 계산 핵심 공식

```python
import numpy as np
from scipy import stats

def simple_allometry(species_bw, species_cl, human_bw=70.0):
    """
    단순 allometry: CL = a × BW^b
    log-log 선형 회귀로 a, b 추정
    """
    log_bw = np.log(species_bw)
    log_cl = np.log(species_cl)
    b, log_a, r, p, se = stats.linregress(log_bw, log_cl)
    a = np.exp(log_a)
    human_cl = a * human_bw ** b
    return {
        "a": a, "b": b, "r2": r**2,
        "human_cl": human_cl,
        "fold_change_from_rat": human_cl / species_cl[species_bw.index(0.25)]
    }

def rule_of_exponents(bw, cl, human_bw=70.0):
    """
    Rule of Exponents (ROE) — Mahmood & Balian 1996
    b < 0.55:  단순 스케일링
    0.55-0.70: 단순 스케일링
    0.71-0.99: MLP 보정 필요
    ≥ 1.00:   뇌중량 보정 필요
    """
    log_bw = np.log(bw)
    log_cl = np.log(cl)
    b, log_a, _, _, _ = stats.linregress(log_bw, log_cl)
    a = np.exp(log_a)
    
    MLP = {"mouse":3.5,"rat":4.5,"rabbit":9,"dog":24,"monkey":30,"human":100}
    BrW = {"mouse":0.41,"rat":1.8,"rabbit":12.1,"dog":118,"monkey":106,"human":1450}
    
    if b < 0.71:
        return a * human_bw**b, "simple", b
    elif 0.71 <= b < 1.00:
        # MLP 보정 (종별 MLP 값이 필요)
        return a * human_bw**b, "MLP_needed", b
    else:
        # BrW 보정
        return a * human_bw**b, "BrW_needed", b

def ivive_clearance(clint_uLmin_mgprot,
                    fu_plasma=0.1,
                    fu_mic=1.0,
                    MPPGL=45,
                    liver_weight_g=1500,
                    QH=70.0):
    """
    IVIVE: In vitro CLint → In vivo CL (Well-stirred liver model)
    
    CLint (µL/min/mg protein) → CL_H (L/h)
    
    공식:
      CLu_int = CLint / fu_mic × MPPGL × liver_weight × 60/1e6
      CL_H = QH × fu × CLu_int / (QH + fu × CLu_int)
    """
    # 단위 변환: µL/min/mg → L/h
    CLu_int = clint_uLmin_mgprot / fu_mic * MPPGL * liver_weight_g * 60 / 1e6
    
    # Well-stirred 모델
    CL_H = (QH * fu_plasma * CLu_int) / (QH + fu_plasma * CLu_int)
    extraction_ratio = CL_H / QH
    
    return {
        "CL_H": CL_H,
        "CLu_int": CLu_int,
        "extraction_ratio": extraction_ratio,
        "high_extraction": extraction_ratio > 0.7
    }
```

## UniPlot 생성 코드

```python
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

def create_uniplot(species_data: dict, predictions: dict, drug: str):
    """
    UniPlot: 종간 PK 비교 Log-Log 산포도 + Allometry 회귀선
    
    species_data: {"mouse": {"bw": 0.02, "CL": 0.05}, "rat": {...}, ...}
    predictions: {"simple": 8.2, "roe": 7.8, "ml_ensemble": 8.5, ...}
    """
    colors = {
        "mouse": "#E85D24", "rat": "#854F0B", "rabbit": "#185FA5",
        "dog": "#0F6E56", "monkey": "#534AB7", "human": "#1a1a18"
    }
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["CL — Log-Log Allometry", "Method Comparison"],
        specs=[[{"type": "scatter"}, {"type": "scatter"}]]
    )
    
    # ── Panel 1: Log-Log Allometry ───────────────────────────────────────
    bw_list, cl_list = [], []
    for species, data in species_data.items():
        if species == "human":
            continue
        fig.add_trace(go.Scatter(
            x=[data["bw"]], y=[data["CL"]],
            mode="markers+text",
            name=species.capitalize(),
            marker=dict(color=colors.get(species,"gray"), size=12),
            text=[species.capitalize()],
            textposition="top center"
        ), row=1, col=1)
        bw_list.append(data["bw"])
        cl_list.append(data["CL"])
    
    # Allometry 회귀선
    if len(bw_list) >= 2:
        log_bw = np.log(bw_list)
        log_cl = np.log(cl_list)
        from scipy import stats
        b, log_a, _, _, _ = stats.linregress(log_bw, log_cl)
        a = np.exp(log_a)
        
        bw_range = np.logspace(-2, 2, 100)  # 0.01 ~ 100 kg
        cl_range = a * bw_range ** b
        fig.add_trace(go.Scatter(
            x=bw_range, y=cl_range,
            mode="lines",
            name=f"Allometry (b={b:.3f})",
            line=dict(color="#534AB7", dash="dash")
        ), row=1, col=1)
        
        # 인간 예측점
        human_pred = a * 70.0 ** b
        fig.add_trace(go.Scatter(
            x=[70.0], y=[human_pred],
            mode="markers",
            name="Human prediction",
            marker=dict(color="#1a1a18", size=16, symbol="star")
        ), row=1, col=1)
    
    # ── Panel 2: 방법별 예측값 비교 ────────────────────────────────────
    methods = list(predictions.keys())
    values  = list(predictions.values())
    
    colors_methods = ["#534AB7","#0F6E56","#854F0B","#185FA5","#E85D24",
                      "#A32D2D","#72243E","#3B6D11","#993C1D","#5F5E5A"]
    
    fig.add_trace(go.Bar(
        x=methods, y=values,
        name="Method predictions",
        marker_color=colors_methods[:len(methods)],
        text=[f"{v:.2f}" for v in values],
        textposition="outside"
    ), row=1, col=2)
    
    # 앙상블 평균선
    ensemble_mean = np.mean(values)
    fig.add_hline(y=ensemble_mean, line_dash="dash", line_color="black",
                  annotation_text=f"Ensemble mean: {ensemble_mean:.2f}",
                  row=1, col=2)
    
    fig.update_xaxes(type="log", title_text="Body Weight (kg)", row=1, col=1)
    fig.update_yaxes(type="log", title_text="CL (L/h)", row=1, col=1)
    fig.update_xaxes(tickangle=45, row=1, col=2)
    fig.update_yaxes(title_text="CL (L/h)", row=1, col=2)
    
    fig.update_layout(
        title=f"UniPlot — {drug.capitalize()} Allometry Ensemble",
        height=500, width=1100, showlegend=True
    )
    
    return fig

# 사용 예시
species_data = {
    "mouse":  {"bw": 0.02, "CL": 0.1},
    "rat":    {"bw": 0.25, "CL": 0.5},
    "rabbit": {"bw": 2.5,  "CL": 2.1},
    "dog":    {"bw": 10.0, "CL": 4.8},
    "monkey": {"bw": 5.0,  "CL": 3.2},
}
predictions = {
    "simple": 7.82, "fixed_b075": 8.12, "roe": 7.65,
    "mlp": 7.90, "brw": 8.30, "two_species": 7.71,
    "fu_corrected": 7.55, "ivive": 8.05, "mpbpk": 7.88, "ml": 8.01
}
fig = create_uniplot(species_data, predictions, "busulfan")
fig.write_html("uniplot_busulfan.html")
```

## 앙상블 선택 기준 (요약)

| 방법 | 적용 조건 | 예측 오차 (fold) |
|------|---------|--------------|
| Simple | 3종+ 데이터, 선형 log-log | ±3배 |
| Fixed b=0.75 | 이론적 근거 우선 | ±4배 |
| ROE | b값 범위로 자동 선택 | ±2배 (선택 시) |
| MLP 보정 | b = 0.71–0.99 | ±2배 |
| fu 보정 | 단백결합 데이터 있을 때 | ±2배 |
| IVIVE | in vitro CLint 있을 때 | ±3배 |
| ML 앙상블 | 4종+ 데이터 | ±2배 (LOO 검증) |
| **Ensemble 최종** | 모든 방법 중위수 | **±2배 목표** |
