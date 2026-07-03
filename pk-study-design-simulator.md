---
name: pk-study-design-simulator
description: >
  약동학(PK) 임상시험 연구설계 시뮬레이터 스킬.
  약물 특성(PK 파라미터, 전임상/선행임상 데이터)을 입력하면
  최적 채혈 설계(sparse/dense), 대상 환자 수, 용량 범위, 투여 일정을 
  Monte Carlo 시뮬레이션으로 탐색하고 최적 설계안 3가지를 비교 제안함.
  
  Interactive HTML 대시보드로 설계 파라미터를 조작하며 실시간으로
  예측 정밀도(RSE%), 파워, 목표 AUC 달성률을 확인할 수 있음.

  다음 요청에 반드시 사용:
  - "PK 연구 설계", "임상시험 설계", "채혈 스케줄 최적화"
  - "몇 명 환자가 필요한가", "sample size for PK study"
  - "sparse sampling 설계", "최적 채혈 시점", "D-optimal design"
  - "전임상 데이터로 임상 설계", "First-in-human PK 설계"
  - "약물 특성 입력해서 연구 설계", "study design simulation"
  - "1상 PK 설계", "2상 PK/PD 설계", "3상 집단약동학 설계"
  - "PopPK 서브스터디 설계", "희소 채혈 최적화"
  - "어떤 용량 범위로 임상할까", "FIH dose escalation 설계"
  - 약물명 + "연구 설계해줘", "임상 디자인 제안해줘" 류의 요청

  출력물:
  1. Interactive HTML 연구설계 대시보드 (슬라이더 + 시뮬레이션)
  2. 최적 설계안 3개 비교표 (채혈 시점, N, 예측 정밀도)
  3. Power calculation 요약
  4. 규제 준수 체크리스트 (FDA/ICH E14, EMA 가이드라인)

compatibility:
  - requires: bash_tool, create_file, present_files
  - optional: web_search (규제 가이드라인 최신화), PubMed (선행 PK 문헌)
---

# PK Study Design Simulator

## 개요

약물의 전임상 또는 임상 선행 데이터를 입력받아 **최적 PK 연구 설계**를 시뮬레이션으로 도출합니다.
핵심 원칙: **D-optimality + Practical constraints + Regulatory compliance**

---

## 에이전트 실행 파이프라인 (8단계)

### Stage 0 — 입력 정보 수집 (Interactive)

```
[필수 입력]
  A. 약물명 / 모달리티 (SM / Protein / ADC / 세포치료제)
  B. 투여 경로 (IV bolus / IV infusion / PO / SC)
  C. 개발 단계 (전임상 → FIH / Ph1 / Ph2 / Ph3 PopPK 서브스터디)
  D. 선행 데이터 유형:
     - 전임상만: PK 파라미터 스케일링 (동물 → 인간)
     - Phase 1 결과 있음: 직접 파라미터 입력
     - 문헌 값만: 유사 약물 기반 추정

[선택 입력 — 있으면 정밀도 향상]
  E. 예상 CL (L/h/70kg) ± 불확실성 범위
  F. 예상 V (L/70kg) ± 불확실성 범위  
  G. 예상 t½ (h) — 또는 CL + V에서 자동 계산
  H. IIV 크기 (CV% for CL, V) — 없으면 DB 기본값 적용
  I. 목표 PK 지표 (AUC 목표, Ctrough 목표, TI 정보)
  J. 실용적 제약 (최대 채혈 횟수/일, 병원 방문 가능 시점)
  K. 특수 집단 (소아, 노인, 신/간기능 저하)

→ 없는 정보는 references/pk-prior-database.md + web_search로 보완
→ 불확실성 클수록 → 더 보수적(큰 N) 설계 권고
```

### Stage 1 — 약물 PK 프로파일 구축

```
입력 경로별 처리:

[전임상 데이터만 있는 경우]
  동물 → 인간 스케일링:
  CL_human = CL_animal × (BW_human/BW_animal)^0.75 × R_CL
  V_human  = V_animal  × (BW_human/BW_animal)^1.00 × R_V
  
  종별 보정 인자 (R_CL):
    Rat→Human: 0.40 (보수적), 평균 0.25
    Dog→Human: 0.65
    Monkey→Human: 0.85
  
  단백결합 보정:
  CL_u_human = CL_total × (fu_human / fu_animal)

[Phase 1 데이터 있는 경우]
  NCA 또는 PopPK 결과 직접 사용
  IIV 추정값 그대로 적용

[문헌/유사 약물 기반]
  references/pk-prior-database.md에서 클래스별 평균 파라미터 적용
  불확실성 ×2~5배로 설정
```

### Stage 2 — 연구 목적 및 설계 유형 분류

```
[FIH / Phase 1 PK]
목적: 안전성 + PK 특성화 + MTD/RP2D 결정
설계 유형: Single-ascending dose (SAD) + Multiple-ascending dose (MAD)
채혈: Dense (골조) — 10~12 시점
코호트당 N: 3+3 또는 accelerated titration
PK 목표: CL, V, t½, Cmax, AUC 정확도 (%RSE <30%)

[Phase 2 PK/PD]
목적: E-R 관계 확립, 용량 최적화
설계 유형: Fixed dose groups + PK 서브스터디
채혈: Sparse (3~5 시점) + optional dense 코호트
N: 30~100명 (PD endpoint 기반)
PK 목표: PopPK 파라미터 %RSE <20%

[Phase 3 PopPK 서브스터디]
목적: 공변량 영향 정량화, 특수 집단 외삽
설계 유형: Opportunistic sparse sampling
채혈: 2~4 시점/환자 (pre-dose, 1~2 Css 시점)
N: 100~300명 (다양한 공변량 분포)
PK 목표: CL, V의 공변량 효과 %RSE <30%
```

### Stage 3 — D-Optimal 채혈 시점 설계

```
D-Optimality 원칙:
  정보 행렬 M(ξ) = Σ F'(tᵢ) × W(tᵢ) × F(tᵢ)
  det[M(ξ)]를 최대화하는 시점 집합 {t₁, t₂, ..., tₙ} 탐색

실용적 근사 알고리즘 (references/d-optimal-algorithm.md):
  1. 1구획 IV 주입: 필수 시점 = 분포 초기, 제거 중기, 제거 말기
  2. 1구획 PO:    필수 시점 = 흡수 상승, Cmax, 제거 2점
  3. 2구획 IV:   필수 시점 = α-phase, β-phase 각 2점

희소 채혈 최적 시점 (Sparse — 4 시점 기준):
  t* ≈ [0.2×t½, 0.7×Tmax, 1.5×t½, 4×t½]
  
  또는 PFIM/PopED 소프트웨어 결과 준용:
  references/pfim-results-library.md 참조

실용적 제약 통합:
  - 병원 방문 제약 (예: 1h, 2h, 4h, 6h, 24h만 가능)
  - 최대 채혈량 제한 (소아: 3 mL/kg, 성인: 500 mL/study)
  - 검체 처리 창 (원심분리 전 최대 30분)
```

### Stage 4 — 샘플 크기 결정 (N 계산)

```
PopPK 기반 N 계산 (simulation approach):

알고리즘:
  FOR N in [20, 30, 50, 75, 100, 150, 200]:
    1. N명 가상 환자 생성 (공변량 분포 반영)
    2. 각 환자에게 선택된 채혈 시점 적용
    3. FOCE 추정 → %RSE 계산 (CL, V 각각)
    4. 목표 %RSE 달성 여부 확인
  최소 달성 N → 권고 N (20% 탈락률 보정: N_enroll = N/0.8)

목표 정밀도 기준:
  주요 파라미터 (CL, V): %RSE ≤ 20%
  IIV (ω²): %RSE ≤ 40%
  잔차 변동 (σ²): %RSE ≤ 30%
  공변량 효과 (θ_cov): %RSE ≤ 30%

파워 계산 (공변량 검출력):
  최소 검출 가능 효과: θ_cov = 0.20 (±20% CL 차이)
  목표 파워: ≥ 80% (α = 0.05, ΔOFV ≥ 3.84 기준)
  → 공변량별 N 요구량 계산
```

### Stage 5 — 용량 범위 및 일정 최적화

```
[FIH Dose Escalation]
  시작 용량 원칙:
    MRSD (Minimum Risk Starting Dose) = NOAEL_animal / (10×HED인자)
    또는 PAD (Pharmacologically Active Dose) 기반
  
  에스컬레이션 규칙:
    전통적 3+3: 용량 제한독성(DLT) 0/3 → 다음 코호트
    加속: EWOC (Escalation With Overdose Control)
    Bayesian: BLRM (Bayesian Logistic Regression Model)
  
  최대 에스컬레이션 단계 (기본값):
    ×2배씩 → 예상 MTD 50% 도달 후 ×1.5배 전환

[반복 투여 일정]
  Steady-state 도달 시점: ~5×t½
  PK 샘플링 코호트:
    - 1일차 (단회): PK 프로파일
    - X일차 (SS): Css,min + Css,max
  
  Q 최적 투여 간격:
    τ_opt = t½ × ln(2) (fluctuation 최소화)
    또는 목표 Ctrough / Cmax 비율 기반 계산
```

### Stage 6 — 시뮬레이션 실행 및 설계안 비교

```
시뮬레이션 결과로 설계안 3개 생성:

[Design A — 최소 부담 (Minimal Burden)]
  채혈 시점: 3~4개 (실용적 최소)
  N: 최소 통계적 요건 충족
  비용: 낮음
  정밀도: 주요 파라미터 %RSE 25~30%
  권장 상황: Phase 3 기회적 sampling, 소아/노인 연구

[Design B — 균형 (Balanced) ★ 권고]
  채혈 시점: 5~7개 (sparse + 핵심 dense)
  N: 목표 정밀도 기준 충족
  비용: 중간
  정밀도: CL/V %RSE ≤ 20%, IIV %RSE ≤ 35%
  권장 상황: Phase 2 PK/PD, 일반적 PopPK

[Design C — 최대 정보 (Maximum Information)]
  채혈 시점: 10~12개 (전통적 dense)
  N: 공변량 효과 검출 최적화
  비용: 높음
  정밀도: CL/V %RSE ≤ 15%, IIV %RSE ≤ 25%
  권장 상황: FIH PK, 복잡한 PK 특성화, NDA 제출용

비교 지표:
  - %RSE (CL, V, IIV, σ²)
  - 목표 AUC 달성률 (% in target range)
  - 총 채혈량 (mL)
  - 환자 부담 지수 (방문 횟수 × 시간)
  - 예상 소요 비용 지수
```

### Stage 7 — 규제 요건 체크리스트

```
references/regulatory-checklist.md 기반 자동 평가:

[FDA 요건]
  □ Guidance for Industry: Population PK (1999)
  □ FDA Model-Informed Drug Development (MIDD) 2017
  □ ICH E4 (용량 반응 연구)
  □ ICH E5 (민족적 요인의 수용성)
  □ ICH E14 (QT/QTc 연장, 약물에 따라)
  □ 소아: FDA PREA/BPCA 준수
  □ 노인: FDA 노인 연구 가이드라인

[EMA 요건]
  □ EMA Population PK Guideline (2007)
  □ EMA DDI Studies Guideline (2012)
  □ Paediatric Investigation Plan (PIP) 해당 시

[분석법 요건]
  □ FDA Bioanalytical Method Validation (2018)
  □ EMA Bioanalytical Method Validation (2011)
  □ 검출한계, 정량한계, 선형 범위 설정
  □ 매트릭스 효과 평가

[데이터 관리]
  □ 21 CFR Part 11 (전자 기록)
  □ CDISC SDTM/ADaM 데이터 형식
  □ NONMEM/Phoenix 소프트웨어 검증
```

### Stage 8 — 출력물 생성

```
[Interactive HTML 대시보드 생성]
  → 다음 요소 포함:
  
  Panel 1: 약물 PK 프로파일 입력 폼
    - CL, V, ka, t½ 슬라이더
    - IIV 크기 조정
    - 투여 경로 선택
  
  Panel 2: 설계 파라미터 조정
    - 채혈 시점 드래그-앤-드롭
    - N 슬라이더
    - 용량 범위 설정
  
  Panel 3: 실시간 시뮬레이션 결과
    - 농도-시간 곡선 (중위수 + 90% PI)
    - %RSE 실시간 업데이트
    - 목표 AUC 달성률 바 차트
  
  Panel 4: 설계안 비교 테이블 (A vs B vs C)
  
  Panel 5: 규제 체크리스트 (자동 평가)

[xlsx 요약 보고서]
  Sheet 1: 입력 파라미터 요약
  Sheet 2: 시뮬레이션 결과 (N별 %RSE)
  Sheet 3: 설계안 3개 비교
  Sheet 4: 채혈 일정표 (달력 형식)
  Sheet 5: 규제 체크리스트

→ 출력 후 present_files로 공유
→ 추가 조정 필요 시 슬라이더 값 변경 요청 안내
```

---

## 개발 단계별 설계 표준 (빠른 참조)

| 단계 | 설계 유형 | 채혈 수/환자 | 최소 N | 목표 정밀도 |
|------|---------|------------|-------|-----------|
| 전임상 PK | 위성 그룹 | 10~12 | 3/성별/용량 | 개인 NCA |
| FIH (SAD) | Dense + 코호트 에스컬레이션 | 10~12 | 6/코호트 | CL/V RSE<30% |
| FIH (MAD) | Dense Day1 + Dense DayN(SS) | 8~10 | 6/코호트 | Css, 축적 |
| Ph1b 용량확장 | Dense 코호트 + sparse 보조 | 6~8 | 12~20 | E-R 초기 |
| Ph2 PK/PD | Sparse(PopPK) + dense 서브그룹 | 4~6 | 30~50 | CL/V RSE<25% |
| Ph3 PopPK | Opportunistic sparse | 2~4 | 100~200 | 공변량 RSE<30% |
| 특수집단 | Dense cross-over 또는 매칭 | 6~10 | 8~12/그룹 | 그룹 간 비교 |

---

## 전임상→임상 스케일링 자동 처리

```
입력: 동물 PK 데이터
출력: 예측 인간 PK + 불확실성 범위

알고리즘:
  1. 단순 allometric: CL_h = a × BW^0.75
  2. Rule-of-Exponents (ROE): t½ 보정 포함
  3. IVIVE (In Vitro – In Vivo): 마이크로솜 CLint + fu,mic → CL_h
  4. PBPK 스케일링: 조직 분배 계수 (Kp) 기반

불확실성 보고:
  예측 CL_h 범위: 중위값 × [1/3 ~ 3] (통상적 오차)
  설계에 보수적 값 (큰 CL, 작은 V) 적용 권고
```

---

## 핵심 참조 파일

```
references/
├── pk-prior-database.md          ← 약물 클래스별 기본 PK 파라미터
├── d-optimal-algorithm.md        ← D-최적 설계 알고리즘 상세
├── pfim-results-library.md       ← PFIM/PopED 시뮬레이션 결과 라이브러리
├── regulatory-checklist.md       ← FDA/EMA 규제 요건 체크리스트
├── sample-size-tables.md         ← 시나리오별 N 계산 참조표
├── allometric-scaling.md         ← 동물→인간 스케일링 공식
├── dose-escalation-rules.md      ← FIH 에스컬레이션 규칙 (3+3, BLRM)
└── html-dashboard-template.md    ← Interactive HTML 출력 템플릿
```

---

## 투여경로별 PK 모델 자동 선택 (Stage 2 확장)

에이전트는 투여경로 입력을 받으면 `references/route-pk-models.md`를 로드하여
해당 경로의 기본 구조 모델, NONMEM 코드, 채혈 전략을 자동 적용합니다.

### 경로 → 모델 매핑 규칙

```
IV bolus       → 1CMT/2CMT, ADVAN1/3, TRANS2/4
IV infusion    → 1CMT/2CMT + RATE 컬럼, Css 계산
IM 즉시방출    → 1구획 1차 흡수, ADVAN2 TRANS2, F≈1
IM depot       → Krel + Ka 순차 모델, ADVAN2 + 0차 방출
SC 즉시방출    → IM 즉시와 동일, F 추정 포함
SC depot/implant → 0차 방출 + 1차 흡수
Biologics SC   → 2구획 + 림프 흡수 경로, 느린 Ka
PO (경구)      → 1차 흡수 1CMT/2CMT, F·Tlag 추정
SL (설하)      → PO 모델 + 매우 큰 Ka, F_SL 추정
Buccal (협측)  → SL과 유사, 0차/1차 복합 가능
Pulmonary/Nasal → 3구획 deposition 모델 (폐+GI+혈중)
Transdermal    → 0차 방출 + 피부 저장층, 제거 후 농도
Vaginal/Rectal → PO 유사 + 1st-pass 부분 회피
Ocular         → 다구획 안구 모델, 국소/전신 이중 추적
```

### 투여경로별 필수 채혈 시점 자동 계산

```python
def optimal_sampling(route, t_half, n_timepoints):
    """
    투여경로와 t½ 기반으로 D-optimal 채혈 시점 반환
    """
    early_dense = {
        'IV_bolus':    [0.08, 0.25, 0.5],   # 분포상 포착
        'IV_infusion': [0.5, 0.8, 1.0],     # Tinf 대비 상대 시점
        'IM_SC':       [0.5, 1.0, 2.0],     # 흡수상
        'Biologics_SC':[12, 24, 48],         # 느린 흡수 (h)
        'PO':          [0.5, 1.0, 1.5],     # ka 추정
        'SL_Buccal':   [0.08, 0.25, 0.5],   # 빠른 흡수
        'Pulmonary':   [0.08, 0.25, 0.5],   # 흡입 후 즉시
        'Transdermal': [4, 8, 12],           # 느린 경피 흡수
        'Vaginal_Rectal':[0.5, 1.0, 2.0],
        'Ocular':      [0.08, 0.5, 1.0],
    }
    late_sparse = [1*t_half, 2*t_half, 4*t_half]
    return early_dense[route] + late_sparse
```
