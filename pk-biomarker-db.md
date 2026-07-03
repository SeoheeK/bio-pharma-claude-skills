---
name: pk-biomarker-db
description: >
  PK/PD 모델 설계 및 신약 개발 시 Target Biomarker 추천 스킬.
  약물 표적(target), 작용기전(MoA), 치료 영역별로 활성/저해가 예측되는
  약동력학적(PD) 바이오마커를 추천하고, PK/PD 모델 구조를 제안함.
  ClinicalTrials.gov 프로토콜 + FDA 심사문서 + 주요 논문 기반.

  다음 요청에 반드시 사용:
  - "이 약물의 PD 마커 추천해줘", "바이오마커 뭘 써야 해?"
  - "target engagement 마커", "약력학 마커 추천"
  - "PK/PD 모델 설계할 때 어떤 마커?", "효능 예측 마커"
  - "독성 바이오마커", "안전성 마커 추천"
  - "면역항암제 바이오마커", "항암제 PD 마커"
  - "VEGF", "IL-6", "ctDNA", "종양 크기" 관련 PK/PD 모델 설계
  - "임상 PD 마커 선택", "phase 1 바이오마커 전략"
  - "target biomarker activator inhibitor 예측"

  출력물:
  1. 약물/표적별 권장 PD 바이오마커 목록 (활성/저해 분류)
  2. PK/PD 모델 구조 제안 (직접/간접 반응 모델)
  3. 마커별 측정 방법 + 채혈 시점
  4. 임상 단계별 바이오마커 전략
  5. NONMEM PK/PD 모델 코드 스니펫

compatibility:
  - requires: create_file, present_files
  - optional: web_search (최신 마커 업데이트), PubMed:search_articles
---

# PK/PD Biomarker DB

## 바이오마커 분류 체계

```
Level 1: 표적 참여 마커 (Target Engagement, TE)
  → 약물이 표적에 결합했는지 직접 측정
  예: PD-1 수용체 점유율 (RO), BCR-ABL 인산화, VEGFR 인산화

Level 2: 약력학 마커 (Pharmacodynamic, PD)
  → 표적 참여 결과로 나타나는 하위 신호 변화
  예: IL-6↓(anti-IL-6R), sVEGF↑(Bevacizumab), ANC↓(항암제)

Level 3: 임상 활성 마커 (Efficacy Biomarker)
  → 효능과 상관있는 종양/질병 반응 지표
  예: ctDNA↓, 종양 크기 (RECIST), PFS, ORR

Level 4: 독성/안전성 마커 (Safety Biomarker)
  → 부작용 예측 또는 감시
  예: Troponin↑(심독성), ALT↑(간독성), ANC↓(골수억제)
```

---

## DOMAIN 1: 면역항암제 (Immune Checkpoint Inhibitors)

### 1-A: 항PD-1/PD-L1 (Nivolumab, Pembrolizumab, Atezolizumab...)

#### Level 1 — Target Engagement 마커
```
[PD-1 수용체 점유율 (RO) — 핵심 TE 마커]
  측정: 말초혈 T세포 (flow cytometry)
        항약물 항체 vs 레퍼런스 항체 경쟁 결합
  PK/PD 관계:
    RO(t) = C(t) / (C(t) + EC50_RO)     ← Emax 모델
    EC50_RO ≈ 0.1–1.5 µg/mL (pembrolizumab)
  임상 의의:
    RO ≥ 80% → 최대 효과 포화 (FDA DOGE 기준 논의)
    용량 선택 근거로 NCT02573259 활용

  NONMEM PK/PD 코드:
    EC50_RO = THETA(X)
    EFFECT  = C/(C + EC50_RO)  ; 0–1 scale
    RO_obs  = EFFECT * 100 + EPS(2)
```

#### Level 2 — Downstream PD 마커
```
[활성화 T세포 마커]
  - CD8+ T세포 증가 ✅ 예측: 치료 반응 예측
  - Ki67+CD8+ T세포 (증식 중인 effector)
  - 혈중 IL-2, IFN-γ ↑ (사이토카인 방출)
  측정: PBMC 분리 + FACS / ELISA
  시점: 기저, 1주기 중간, 2주기 전

[면역 억제 마커]
  - sCD163 (M2 대식세포 활성) ↓ 예측: 면역억제 해제
  - TGF-β ↓ (면역억제 사이토카인)
  - IL-10 ↓

[조직 마커 (종양 미세환경)]
  - 종양 침윤 림프구 (TIL) density ↑ → 반응 예측
  - PD-L1 TPS/CPS ≥ 50% → 반응률 ↑ (but 완전하지 않음)
  - CD8:Treg 비율 ↑
  - TMB (Tumor Mutational Burden) ≥ 10 mut/Mb → MSI-H 연관
```

#### Level 3 — 효능 마커
```
[ctDNA (circulating Tumor DNA)] ← 현재 가장 유망
  - 치료 후 4주 ctDNA ↓ ≥ 50%: PFS 연장 예측
  - ctDNA 증가: 종양 재발/내성 조기 감지
  - ddPCR 또는 Next-Gen Sequencing
  모델: ctDNA ~ f(종양 크기 × 이탈률)

[종양 크기 (Tumor Growth Inhibition Model, TGI)]
  TGI PK/PD 모델:
    dTS/dt = kg×TS - ks×C×TS    (Simeoni 모델)
    kg: 종양 성장속도
    ks: 약물-종양 억제 상수
    C: 혈중 약물 농도 (PK 연동)

[sVEGF (Soluble VEGF)] — Bevacizumab 병용 시
  sVEGF ↑ 의역: anti-VEGF로 인한 VEGF 포획
```

#### Level 4 — 독성 마커
```
[면역 관련 부작용 (irAE) 예측]
  - IL-17 ↑ → 피부 irAE 위험
  - ANA, anti-dsDNA → 자가면역 위험
  - 기저 Thyroid Ab → 갑상샘 irAE 위험
  - Eosinophil ↑ → 중증 irAE 위험 신호

[사이토카인 방출 증후군 (CRS) — CAR-T]
  - IL-6 ↑↑↑ → CRS 중증도 예측 (Tocilizumab 치료 기준)
  - CRP, Ferritin ↑
  - Grade ≥ 3 CRS 예측 모델: IL-6 > 1000 pg/mL
```

---

## DOMAIN 2: 표적항암제 (Targeted Therapy)

### 2-A: VEGF/VEGFR 표적제 (Bevacizumab, Sunitinib, Sorafenib...)

#### PD 마커
```
[활성(증가) 예측 마커 ← 표적 결합으로 증가]
  sVEGF (혈청 유리 VEGF) ↑↑
    → Bevacizumab: 2–4배 증가 (Fc 결합 → 포획)
    → 포획된 VEGF 복합체 → CL 감소 → TMDD 모델 필요
    
  SDF-1α (CXCL12) ↑
    → 항혈관신생 치료 후 골수 유래 세포 보상 반응
    
  PlGF (Placenta Growth Factor) ↑
    → VEGF 억제 후 보상적 증가 → 내성 신호

[저해(감소) 예측 마커 ← 약물 효과로 감소]
  sKDR (soluble VEGFR2) ↓
    → VEGFR2 인산화 억제 → Sunitinib 효과 지표
    
  혈관 내피 세포 증식 (PECAM-1/CD31) ↓
    → 종양 혈관 신생 억제

  혈압 ↑ (독성이지만 PD 마커)
    → VEGFR 억제 → NO 감소 → 혈압 상승
    → Sunitinib 효능 대리 마커!

PK/PD 모델:
  TMDD (Target-Mediated Drug Disposition):
    dRtot/dt = ksyn - kdeg×Rtot - kon×C×Rtot + koff×RC
    dRC/dt   = kon×C×Rtot - (koff+kint)×RC
    C_obs    = (Ctot - RC) / V
```

### 2-B: BCR-ABL 표적제 (Imatinib, Dasatinib, Nilotinib...)

#### PD 마커
```
[저해 마커 (표적 저해)]
  BCR-ABL 키나제 활성 ↓
    → pCRKL (phospho-CRKL) ↓: 매우 예민한 TE 마커
    → pSTAT5 ↓
    → Western blot 또는 flow cytometry

  CML 반응 지표 (효능):
    BCR-ABL1 IS % ↓ (국제표준, 분자생물학적 반응)
    → MR4.0 (BCR-ABL < 0.01%) = 깊은 분자 반응

  ANC (Absolute Neutrophil Count) ↓
    → 골수억제 — 독성이지만 PD 마커
    PK/PD 모델: Semi-mechanistic myelosuppression
      (Friberg et al. J Clin Oncol 2002)

[NONMEM Myelosuppression 코드]
  $PK
    MTT  = THETA(1)    ; Mean transit time
    SLOPE= THETA(2)    ; Drug effect slope
    CIRC0= THETA(3)    ; Baseline ANC

  $DES
    EDRUG  = SLOPE*C
    DADT(1)= CIRC0*(1-EDRUG)/MTT - A(1)/MTT  ; Prol
    DADT(2)= A(1)/MTT - A(2)/MTT              ; Trans1
    DADT(3)= A(2)/MTT - A(3)/MTT              ; Trans2
    DADT(4)= A(3)/MTT - A(4)/MTT              ; Trans3 → Circ
    CIRC   = A(4)
```

### 2-C: EGFR 표적제 (Erlotinib, Gefitinib, Osimertinib...)

#### PD 마커
```
[저해 마커]
  pEGFR ↓ (종양 생검/세포주)
  pERK ↓, pAKT ↓ (하위 신호)
  HER2, HER3 발현 변화

[저항성 마커 — 2차/3차 돌연변이]
  T790M 돌연변이 ↑ → 1세대 EGFR TKI 내성 (Osimertinib 적응)
  C797S 돌연변이 ↑ → 3세대 내성
  ctDNA에서 돌연변이 동적 추적

[피부 부작용 — PD 마커]
  Acneiform rash ← EGFR 억제 → 피부 부작용
  역설적 효능 예측: 피부독성 경험자 = 더 나은 효능?

[효능 예측]
  EGFR 돌연변이 (Ex19del, L858R) → 반응 예측 마커
  FISH amplificate → 고반응 예측
```

---

## DOMAIN 3: 항염증/자가면역제 (Anti-inflammatory)

### 3-A: Anti-IL-6/IL-6R (Tocilizumab, Sarilumab...)

#### PD 마커
```
[저해 마커 (IL-6 신호 억제)]
  CRP (C-Reactive Protein) ↓↓↓  ← 가장 빠른 PD 마커 (24–48h)
    → Tocilizumab 후 CRP 정상화: 효능 확인
    → PK/PD 모델: Indirect Response Model
    
  SAA (Serum Amyloid A) ↓
  ESR (Erythrocyte Sedimentation Rate) ↓
  Fibrinogen ↓
  IL-6 혈중 농도 ↑ (역설적)
    → IL-6R 차단 → IL-6 반감기 증가 → 혈중 축적
    → 총 IL-6 ↑지만 활성 IL-6 신호 ↓

[활성화 마커 (저해 해제 = 증가)]
  sIL-6R ↑ (가용형 IL-6 수용체 증가)
    → 수용체 내재화 억제로 shed

[NONMEM Indirect Response PD 모델]
  dCRP/dt = kin - kout × CRP
  kin0    = kout × CRPbaseline  ; 기저 상태
  kin(t)  = kin0 × (1 - Imax × C^n / (IC50^n + C^n))
  
  ; Imax: 최대 억제율 (~1.0)
  ; IC50:  반최대 억제 농도
  ; n:     Hill 계수 (통상 1~2)
```

### 3-B: Anti-TNF (Adalimumab, Infliximab...)

#### PD 마커
```
[저해 마커]
  TNF-α ↓ (유리 TNF)
    → 항TNF 제제 결합 → 유리 TNF 감소
    → BUT 총 TNF+약물 복합체 측정 어려움
    
  CRP ↓ (IL-6 이차 감소)
  DAS28 / CDAI / MAYO score 개선 (RA, IBD)
  
  sCD163 ↓ (M2 대식세포 활성 저하)
  IL-1β ↓, IL-17A ↓ (관절 세포)

[내성/면역원성 마커]
  ADA (Anti-Drug Antibody) ↑ → 약물 CL ↑ → 효능 상실
    → ADA 모니터링 필수
    → ADA 형성률: Adalimumab 15–20%, Infliximab 30–60%
  
[Trough 농도 — TDM PD 마커]
  Infliximab Ctrough ≥ 3–7 µg/mL → 점막 치유
  Adalimumab Ctrough ≥ 5–8 µg/mL → 임상 관해

PK/PD 목표 농도:
  약물 Ctrough vs 반응률 E-R 모델:
    Effect(Ctrough) = Emax×Ctrough / (EC50 + Ctrough)
    Emax ≈ 0.85 (최대 85% 반응률)
    EC50 ≈ 2.0 µg/mL (Adalimumab 기준)
```

---

## DOMAIN 4: 항감염제 (Antimicrobials)

### 4-A: β-락탐 항생제 (Penicillin, Cephalosporin, Carbapenem)

#### PD 마커
```
[항균 PD 지표 — Time-Dependent]
  fT>MIC (free drug T above MIC) ← 핵심 PD 인덱스
    → Penicillin/Cephalosporin: fT>MIC ≥ 40–70% τ
    → Carbapenem: fT>MIC ≥ 40% τ
    
  MIC (Minimum Inhibitory Concentration) 측정
  세균 계수 (CFU) 동적 변화 → time-kill curve

PK/PD 통합 모델:
  dB/dt = kg×B - kmax×(C/MIC)^γ / (C50/MIC)^γ + (C/MIC)^γ) × B
  또는 Hill 모델 직접 적용:
  Effect(C) = Emax × C^γ / (EC50^γ + C^γ)

  Pharmacodynamic Target Attainment (PTA):
    P(fT>MIC ≥ 40%) by Monte Carlo Simulation
    → 용량 최적화에 직접 활용

독성 마커:
  SCr ↑ (Imipenem 신독성, 고용량)
  Na+ 장애 (고용량 Penicillin)
```

### 4-B: 글리코펩타이드 (Vancomycin)

#### PD 마커
```
PD 인덱스: AUC24/MIC ≥ 400–600 (2020 FDA 개정 기준)
  → 과거 Ctrough-based TDM에서 AUC-based로 전환

PK/PD 통합:
  AUC_24h = Daily dose / CL_individual (NONMEM 추정)
  PTA at MIC = 1 mg/L: P(AUC ≥ 400) 목표

독성 PD 마커:
  AKI (Acute Kidney Injury):
    SCr ↑ ≥ 0.5 mg/dL or 25% 증가 → 독성 기준
    NGAL ↑ (조기 신독성 마커, SCr보다 빠름)
    Ctrough > 15 mg/L → AKI 위험 ↑
```

---

## DOMAIN 5: 알킬화제 (Busulfan/Cyclophosphamide — HSCT)

#### PD 마커
```
[효능 마커]
  AUC 목표 달성 여부 (TDM PK/PD)
  생착 (Engraftment) 속도 → Day+14, +28 ANC
  혼합 키메리즘 (Mixed Chimerism) → Day+30, +100 (STR 분석)
  MRD (Minimal Residual Disease) → 분자 반응

[독성 마커]
  VOD/SOS (정맥 폐쇄성 질환):
    - Bilirubin ↑ + 체중 증가 + 간종대 → Baltimore 기준
    - PAI-1 (Plasminogen Activator Inhibitor-1) ↑
    - ADAMTS13 ↓ (초기 VOD 마커)
    - Busulfan AUC_total > 20,000 µmol·min/L → 위험 증가

  신경독성:
    - 발작 위험 → Phenytoin 예방 시 CL↑ → AUC↓ → 생착 실패
    - Busulfan Ctrough > 1.5 µmol/L → 신경독성

  폐독성 (IPS):
    - CT 폐 침윤 + KL-6 ↑
    - IL-1Ra, TNF-α ↑ (전처치 관련 사이토카인)

[PK/PD 모델]
  AUC → P(VOD) 로지스틱 모델:
    logit[P(VOD)] = β0 + β1 × AUC_total
    P(VOD) = 1/(1 + exp(-logit))
    
  AUC vs 생착 E-R:
    P(engraft) = Emax × AUC / (EC50 + AUC)
    EC50_engraft ≈ 900 µmol·min/L (단회)
```

---

## DOMAIN 6: 생물의약품 (Biologics — General)

### 항약물항체 (Anti-Drug Antibody, ADA) — 모든 바이오로직스 공통

```
[ADA PK/PD 영향]
  ADA 형성 → CL ↑ → Ctrough ↓ → 효능 상실
  ADA PK 모델:
    CL(t) = CL_base × (1 + θ_ADA × ADA(t))
    또는 혼합효과 모델: ADA 형성자 vs 비형성자 이분

[ADA 모니터링]
  측정 시점: 기저, 매 주기, 의심 시점
  분류: 중화항체(NAb) vs 비중화항체(non-NAb)
  임계값: ADA ≥ LLOQ → 약물 농도 저하 판단

[면역원성 PD 마커]
  - Infliximab ADA(+): Ctrough < 1 µg/mL → 2차 실패
  - 재치료 후 항체 재형성 속도
```

---

## PK/PD 모델 유형 선택 가이드

### 직접 반응 모델 (Direct Response)
```
E(t) = Emax × C(t)^γ / (EC50^γ + C(t)^γ)

적합 상황:
  - 빠른 평형: 수용체 결합/해리 빠름
  - 효과가 즉시 농도와 비례
  예: 항생제 살균 효과 (bactericidal), 국소 마취

NONMEM:
  EMAX = THETA(X)
  EC50 = THETA(Y)
  HILL = THETA(Z)
  E    = EMAX*C**HILL/(EC50**HILL+C**HILL)
  Y    = E + EPS(1)
```

### 간접 반응 모델 (Indirect Response, IDR)
```
모델 유형:
  IDR-1: kin 억제: dR/dt = kin×(1-Imax×C/(IC50+C)) - kout×R
  IDR-2: kout 억제: dR/dt = kin - kout×(1-Imax×C/(IC50+C))×R
  IDR-3: kin 촉진:  dR/dt = kin×(1+Emax×C/(EC50+C)) - kout×R
  IDR-4: kout 촉진: dR/dt = kin - kout×(1+Emax×C/(EC50+C))×R

적합 상황:
  - 효과가 생체 과정을 통해 발현 (단백 합성/분해)
  - 효과-농도 시간 불일치 (hysteresis)
  예: CRP↓(anti-IL6), ANC↓(항암제), INR↑(항응고제)

IDR-1 NONMEM:
  $PK
    KIN  = THETA(1)
    KOUT = THETA(2)
    IMAX = THETA(3)
    IC50 = THETA(4)
    R0   = KIN/KOUT  ; 기저 상태

  $DES
    DADT(1) = KIN*(1-IMAX*C/(IC50+C)) - KOUT*A(1)
    ; A(1) = 반응변수 R (예: CRP)
```

### 전환 구획 모델 (Transit Compartment — 골수억제)
```
Friberg 모델 (골수억제 표준):
  dProl/dt  = ktr×Prol×(1-E(C))×(Circ0/Circ)^γ - ktr×Prol
  dTrans1/dt = ktr×Prol - ktr×Trans1
  dTrans2/dt = ktr×Trans1 - ktr×Trans2
  dTrans3/dt = ktr×Trans2 - ktr×Trans3
  dCirc/dt  = ktr×Trans3 - ktr×Circ

  E(C) = Slope × C  (선형) 또는 Emax 모델
  ktr  = 4/MTT (Mean Transit Time ≈ 5.35일 중성구)

NONMEM: ADVAN6 ODE + $DES 블록
```

### TMDD 모델 (Target-Mediated Drug Disposition)
```
mAb가 표적(수용체)과 결합 → 분포/소실 비선형화

준정상상태 근사 (QSS approximation):
  dL/dt  = kin - (kdeg + kon×(Rtot-LR))×L + koff×LR
  dLR/dt = kon×(Rtot-LR)×L - (koff+kint)×LR
  dRtot/dt = ksyn - kdeg×Rtot - kint×LR

Rtot = 총 수용체, LR = 리간드-수용체 복합체
적합: Bevacizumab/sVEGF, Tocilizumab/sIL-6R, mAb 일반

NONMEM: ADVAN6 + $DES
  KSYN = THETA(X)     ; 수용체 합성속도
  KDEG = THETA(Y)     ; 수용체 분해속도
  KINT = THETA(Z)     ; 복합체 내재화속도
```

---

## 바이오마커 DB 요약 인덱스

| 표적/MoA | Level1 TE 마커 | Level2 PD 마커 (활성↑) | Level2 PD 마커 (저해↓) | Level4 독성 마커 | 추천 PD 모델 |
|--------|-------------|-------------------|-------------------|-------------|-----------|
| PD-1/PD-L1 차단 | PD-1 RO (flow) | CD8+ T세포, IFN-γ, IL-2 | TGF-β, IL-10 | irAE 사이토카인 | Emax (RO), TGI |
| VEGF/VEGFR 억제 | sVEGF↑, sKDR↓ | sVEGF (포획↑), PlGF↑ | pVEGFR2↓, 혈관신생 | 혈압↑, 단백뇨 | TMDD, Indirect |
| BCR-ABL 억제 | pCRKL↓ | BCR-ABL IS %↓, CMR | ANC↓ (골수억제) | 혈소판↓, ANC↓ | Transit (Friberg) |
| EGFR 억제 | pEGFR↓ | ctDNA↓ | pERK↓, pAKT↓ | 피부발진, 설사 | Emax (신호) |
| IL-6R 차단 | sIL-6R↑ | IL-6 (혈중↑) | CRP↓, SAA↓ | 감염 위험 | IDR-1 (CRP) |
| TNF-α 차단 | TNF↓ | CRP↓, DAS28↓ | IL-1β↓ | ADA, 감염 | IDR |
| BCL-2 억제 | pBCL-2↓ | CLL 반응 (절대 림프구수↓) | ANC↓ | 종양용해증후군 | Emax + 안전성 |
| 알킬화제 (BU) | AUC (PK=PD) | 생착↑ (ANC Day+14) | MRD↓ | VOD (빌리루빈↑) | Logistic (VOD P) |
| 항생제 β-락탐 | fT>MIC | 세균 CFU 감소 | 균혈증 지속 | SCr↑ | TKD (Time-kill) |
| 항생제 Vanco | AUC/MIC | MRSA 소실 | AUC/MIC < 400 | NGAL↑, SCr↑ | AUC PD target |
| Busulfan TDM | AUC (PK 직접) | 생착률 | VOD 위험 | 발작, VOD | AUC–반응 E-R |

---

## 바이오마커 측정 시점 권고 (Phase 1 기준)

```
[Phase 1 표준 바이오마커 수집 타임라인]

Day 1 (기저):
  ✅ PBMC (PD-1 RO, lymphocyte subsets)
  ✅ 혈청 사이토카인 패널 (IL-2, IL-6, IL-8, IFN-γ, TNF-α)
  ✅ sVEGF, sVEGFR2 (VEGF 표적 시)
  ✅ CRP, Ferritin, fibrinogen
  ✅ ctDNA (기저 종양 부담)
  ✅ 전혈구검사 (CBC with diff)
  ✅ 간기능, 신기능

Day 8 (Cycle 1, Week 2):
  ✅ PBMC + RO (최대 점유율 확인)
  ✅ 활성화 T세포 (Ki67+CD8+)
  ✅ CRP, 사이토카인

Day 15 (Cycle 1, Week 3):
  ✅ Ctrough (PK 확인 시점)
  ✅ 사이토카인 패널
  ✅ CBC (골수억제 확인)

Day 22 (Cycle 2 전):
  ✅ ctDNA (조기 반응 평가)
  ✅ 안전성 마커 전체

Week 6–8 (Cycle 2–3):
  ✅ 영상 (RECIST) + ctDNA
  ✅ ADA (면역원성 모니터링)
  ✅ 완전 PD 패널

주의: PBMC 처리 2h 이내, 혈청 분리 즉시, -80°C 보관
```
