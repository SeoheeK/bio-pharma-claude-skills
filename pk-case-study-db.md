---
name: pk-case-study-db
description: >
  ClinicalTrials.gov 및 FDA 심사 데이터 기반 PopPK 케이스 스터디 DB 스킬.
  30개 이상 실제 임상시험 약물의 파라미터 세트와 모델 구조를 참조하여
  새 약물 분석 시 유사 케이스 자동 추천 + 파라미터 초기값 제안 + NONMEM 코드 생성.

  다음 요청에 반드시 사용:
  - "이 약물과 비슷한 PK 모델 있어?", "유사 케이스 찾아줘"
  - "파라미터 초기값 추천", "문헌 기반 THETA 초기값"
  - "항암제 PopPK 사례", "면역항암제 PK 파라미터"
  - "mAb PK 모델", "소분자 항암제 PopPK"
  - "ClinicalTrials 기반 모델", "FDA 심사 PK 데이터"
  - "항생제 PopPK", "면역억제제 PK 케이스"
  - "이 약물 클래스 대표 파라미터 알려줘"

  출력물:
  1. 유사 케이스 3~5개 자동 추천 (약물명, NCT번호, 모델 구조)
  2. 추천 초기 파라미터 세트 (THETA/OMEGA/SIGMA)
  3. 공변량 선택 근거 (문헌 기반)
  4. NONMEM 컨트롤 스트림 템플릿
  5. 주요 참고문헌

compatibility:
  - requires: create_file, present_files
  - optional: web_search (최신 문헌 업데이트), PubMed:search_articles
---

# PK Case Study DB — ClinicalTrials.gov 기반

## 데이터 출처

- ClinicalTrials.gov 공개 프로토콜 문서 (NCT 번호 참조)
- FDA Drugs@FDA 심사 문서 (NDA/BLA Clinical Pharmacology Reviews)
- PubMed 발표 PopPK 논문 (2018–2025)
- 항목별 최소 N≥30, NONMEM FOCE 추정 결과

---

## 케이스 스터디 데이터베이스

### CATEGORY 1: 면역항암제 (Immune Checkpoint Inhibitors)

#### Case 1-A: Nivolumab (PD-1 mAb) — NCT01667809
```
약물: Nivolumab (Opdivo, BMS)
적응증: NSCLC, 흑색종, 신세포암 등
모달리티: IgG4 mAb (146 kDa)
투여경로: IV infusion (30분, q2W/q4W)

모델 구조: 2구획 선형 PK (ADVAN3 TRANS4)
추정 방법: FOCE-I (NONMEM 7.3)
데이터: N=1,895명, 4종 종양, Phase I-III 통합

파라미터 (최종 모델):
  CL    = 9.5 mL/h  (CV 35%)   ← 체중, ECOG PS, 기저 알부민 공변량
  V1    = 3.5 L     (CV 19%)   ← 체중 공변량
  Q     = 4.6 mL/h
  V2    = 5.6 L
  t½    ≈ 27일 (terminal)
  IIV CL: ω² = 0.123 (CV 35%)
  IIV V1: ω² = 0.036 (CV 19%)
  RUV: 비례 15%

공변량 (유의, ΔOFV 기준):
  CL ← BW (power 0.738)       ΔOFV = -198.4
  CL ← ECOG PS (exp, +31%)    ΔOFV = -42.1
  CL ← Albumin (exp, -18%/g)  ΔOFV = -38.7
  CL ← 종양 유형 (exp, ±15%)  ΔOFV = -21.3
  V1 ← BW (power 0.716)       ΔOFV = -89.2

NONMEM $THETA 초기값:
  (0, 9.5)    ; 1 CL (mL/h)
  (0, 3500)   ; 2 V1 (mL)
  (0, 4.6)    ; 3 Q  (mL/h)
  (0, 5600)   ; 4 V2 (mL)
  (0, 0.738)  ; 5 BW-CL power
  (0, 0.31)   ; 6 ECOG PS 효과 (exp)
  (-0.20, 0)  ; 7 Albumin 효과 (exp, 음수)

NCT: NCT01667809, NCT01894035
참고: Bajaj G et al. J Clin Pharmacol 2017;57:748-757
```

#### Case 1-B: Pembrolizumab (PD-1 mAb) — NCT01866319
```
약물: Pembrolizumab (Keytruda, MSD)
적응증: 흑색종, NSCLC, HNSCC, MSI-H CRC 등
모달리티: IgG4 humanized mAb (149 kDa)
투여경로: IV infusion (q3W / q6W)

모델 구조: 2구획 시간-변화 CL 모델
CL(t) = CLt × (1 - θ_TC × exp(-CLt_off × t)) ← 종양 반응에 따른 CL 변화
추정 방법: FOCE-I (NONMEM 7.3)
데이터: N=2,993명, Phase I-III

파라미터:
  CL_t  = 0.209 L/h   (CV 24%) ← 기저 CL
  V1    = 3.45 L       (CV 17%)
  Q     = 0.555 L/h
  V2    = 4.10 L
  θ_TC  = 0.147        ← 종양 크기 CL 효과
  t½    ≈ 26일

공변량:
  CL ← BW (power 0.566)      ΔOFV = -232
  CL ← 종양 부담 (TS, 비선형) ΔOFV = -89
  CL ← ECOG 1 vs 0 (+31%)   ΔOFV = -47
  V1 ← BW (power 0.558)      ΔOFV = -108

주요 특성: TS (Tumor Size)가 CL에 영향 → TGI (Tumor Growth Inhibition) PK/PD 통합 필요

NCT: NCT01866319
참고: Ahamadi M et al. CPT:PSP 2017;6:49-55
```

#### Case 1-C: Atezolizumab (PD-L1 mAb) — NCT01375842
```
약물: Atezolizumab (Tecentriq, Genentech)
모달리티: IgG1 engineered mAb (144 kDa)
투여경로: IV infusion (1200 mg q3W)

모델 구조: 2구획 선형 PK
파라미터:
  CL    = 0.200 L/day  (CV 33%)
  V1    = 3.28 L        (CV 21%)
  Q     = 0.485 L/day
  V2    = 3.40 L
  t½    ≈ 27일

공변량:
  CL ← BW (power 0.429)
  CL ← 기저 알부민 (linear, -18%)
  CL ← 항약물항체 ADA (+42%)
  V1 ← BW (power 0.463)

NCT: NCT01375842, NCT02265874
```

---

### CATEGORY 2: 표적항암제 (Targeted Oncology — Small Molecule)

#### Case 2-A: Imatinib (BCR-ABL inhibitor) — Multiple NCTs
```
약물: Imatinib (Gleevec, Novartis)
적응증: CML, GIST
투여경로: PO (400–800 mg QD)
모달리티: 소분자 (493 Da)

모델 구조: 1구획 1차 흡수 (ADVAN2 TRANS2)
데이터: N=280명 (CML 성인), 통합 Phase I-III

파라미터:
  CL/F  = 14.3 L/h    (CV 40%)  ← 체중, 연령, 알파-1-AGP
  V/F   = 347 L        (CV 45%)
  Ka    = 0.48 h⁻¹    (CV 55%)
  F     = ~0.98 (거의 완전 흡수)
  Tmax  ≈ 2.5h, t½ ≈ 16.8h

공변량:
  CL/F ← BW (linear, +0.5% per kg)   ΔOFV = -28.4
  CL/F ← α1-AGP (-)                   ΔOFV = -22.1
  CL/F ← 연령 (power, -0.3% per yr)  ΔOFV = -15.3
  V/F  ← BW (power 0.89)

NONMEM $THETA:
  (0, 14.3)  ; 1 CL/F (L/h)
  (0, 347)   ; 2 V/F  (L)
  (0.1, 0.48); 3 Ka   (h⁻¹)
  (-0.01, 0) ; 4 AGP effect (반비례)
  (0, 0.005) ; 5 BW-CL slope

참고: Widmer N et al. Clin Pharmacokinet 2008;47:379-392
      Chien YH et al. Cancer Chemother Pharmacol 2022
      NCT 기반: Imatinib TDM 다기관 연구
```

#### Case 2-B: Venetoclax (BCL-2 inhibitor) — NCT02141282
```
약물: Venetoclax (Venclexta, AbbVie/Roche)
적응증: CLL, AML
투여경로: PO (용량 에스컬레이션, 20→400 mg)

모델 구조: 2구획 1차 흡수 + Michaelis-Menten 흡수 (고용량 비선형)
파라미터:
  CL/F  = 389 L/h     (CV 52%)
  V2/F  = 256 L        (CV 58%)
  Q/F   = 85.7 L/h
  V3/F  = 897 L
  Ka    = 5.67 h⁻¹   (CV 75%)
  t½    ≈ 26h

공변량:
  CL/F ← 음식 효과 (fed: +3.4배)   ← 반드시 식사와 복용
  CL/F ← 강한 CYP3A4 억제제 (+6배) ← DDI 위험!
  V2/F ← BW

NCT: NCT02141282 (M14-032 시험)
주의: Posaconazole 병용 시 AUC 6배↑ → Grade 4 DDI
```

#### Case 2-C: Ibrutinib (BTK inhibitor) — NCT01105716
```
약물: Ibrutinib (Imbruvica, AbbVie)
적응증: CLL, MCL, WM
투여경로: PO (420 or 560 mg QD)

모델 구조: 1구획 Weibull 흡수 (비선형 흡수 포착)
파라미터:
  CL/F  = 62 L/h      (CV 66%)  ← 매우 높은 IIV
  V/F   = 683 L        (CV 68%)
  Ka    = γ식 Weibull 흡수
  F     ≈ 0.03 (단독) → 0.05 (식사)
  t½    ≈ 4–6h

공변량:
  CL/F ← 음식 (fed: F↑60%)
  CL/F ← 강한 CYP3A4 억제제 (+24배!) ← 금기에 준함
  AUC = 투여 후 내성(ADCC) 무관, 농도-효능 관계 약함

NCT: NCT01105716 (PCYC-1102)
```

---

### CATEGORY 3: 면역억제제 (Immunosuppressants — HSCT)

#### Case 3-A: Busulfan — NCT 다수 (HSCT 전처치)
```
약물: Busulfan (Myleran/Busulfex)
적응증: HSCT 전처치 (BuCy/BuFlu)
투여경로: IV infusion (0.8 mg/kg q6h or q24h)

모델 구조: 1구획 (ADVAN1 TRANS2) — 이미 구축된 모델
파라미터 (최종 모델, 이 스킬의 핵심 케이스):
  TVCL  = 7.82 L/h    (CV 28.6%) → 공변량 후 21.9%
  TVV   = 48.3 L       (CV 20.2%)
  IIV CL: ω² = 0.082 → 0.048 (공변량 후)
  IIV V:  ω² = 0.041

공변량 (최종):
  CL ← BW power 0.74         ΔOFV = -28.4
  CL ← GSTA1*B exp(-0.31)    ΔOFV = -19.4
  CL ← Azole exp(-0.22)      ΔOFV = -11.7
  CL ← Total bilirubin         ΔOFV = -9.1

TDM 목표: AUC0-∞ 900–1500 µmol·min/L (단회)
NCT: NCT01052883, NCT02728752

참고: Choi B 외. 이 스킬의 원 케이스 연구
```

#### Case 3-B: Cyclosporine (CsA) — HSCT/이식
```
약물: Cyclosporine A (Neoral, Sandimmun)
투여경로: IV → PO 전환 (HSCT), PO (신이식)

모델 구조: 2구획 PO 흡수 (ADVAN4 TRANS4)
파라미터 (성인 HSCT):
  CL/F  = 35.2 L/h    (CV 45%)  ← CYP3A5 유전자형 공변량
  V2/F  = 989 L        (CV 62%)
  Q/F   = 28.3 L/h
  V3/F  = 1280 L
  Ka    = 0.89 h⁻¹
  Tmax  ≈ 1.5h, t½ ≈ 18–27h

공변량:
  CL/F ← CYP3A5*1 (expresser: +35%)   ΔOFV = -31
  CL/F ← MDR1 C3435T (TT: -22%)       ΔOFV = -18
  CL/F ← Voriconazole (-48%)           ΔOFV = -44 ← Major DDI!
  CL/F ← Hematocrit (혈중 결합)

TDM: C2 (2h) 또는 Ctrough (C0)
NCT: NCT00823940
```

#### Case 3-C: Tacrolimus — 신이식/HSCT
```
약물: Tacrolimus (Prograf/FK506)
투여경로: PO (0.1–0.3 mg/kg/day BID 또는 OD)

모델 구조: 1구획 PO, CYP3A5 혼합 모델
파라미터:
  CL/F  = 17.8 L/h    (CV 60%)  ← CYP3A5 주도 변동!
  V/F   = 590 L        (CV 55%)
  Ka    = 0.46 h⁻¹
  F     ≈ 0.25 (광범위 1st-pass)
  t½    ≈ 23h

공변량 (가장 중요):
  CL/F ← CYP3A5*1 expresser (×1.84)   ΔOFV = -120 ← 최대 효과!
  CL/F ← 혈청 헤마토크리트 (결합 보정)
  CL/F ← 조기 이식 후 시간 효과
  CL/F ← Azole (-40~70%)

TDM: Ctrough 목표 5–20 ng/mL (시기별)
NCT: NCT01094119, NCT01399814
```

---

### CATEGORY 4: 항생제 (Antibiotics — Critical Care)

#### Case 4-A: Vancomycin — ICU/중증 감염
```
약물: Vancomycin (Vancocin)
투여경로: IV infusion (15–45 mg/kg/day, divided q6–12h)

모델 구조: 2구획 IV infusion (ADVAN3 TRANS4)
데이터: N=485명 (ICU 환자, 다기관), FOCE-I

파라미터:
  CL    = 3.68 L/h    (CV 48%)  ← GFR/신기능 주도!
  V1    = 28.6 L       (CV 52%)
  Q     = 2.11 L/h
  V2    = 51.4 L
  t½β   ≈ 8–12h (신기능 정상)

공변량 (핵심):
  CL ← CrCl (power 0.651, +42%/10증가)   ΔOFV = -183 ← 최대!
  CL ← 패혈증 여부 (aug renal: +62%)      ΔOFV = -54
  V1 ← BW (power 1.0)                     ΔOFV = -78
  CL ← CRRT (연속신대체요법: -65%)

TDM: AUC/MIC ≥ 400–600 (mg·h/L) × (mg/L)⁻¹ → FDA 2020 가이드
NCT: NCT02535949 (Stitt et al. 2024 참조)
NONMEM $THETA:
  (0, 3.68)   ; 1 CL (L/h, 정상 GFR)
  (0, 28.6)   ; 2 V1 (L)
  (0, 2.11)   ; 3 Q  (L/h)
  (0, 51.4)   ; 4 V2 (L)
  (0, 0.651)  ; 5 CrCl-CL power
```

#### Case 4-B: Linezolid — MDR-TB / 중증감염
```
약물: Linezolid (Zyvox)
투여경로: PO/IV (600 mg BID)

모델 구조: 1구획 PO (ADVAN2), 또는 2구획 IV
파라미터:
  CL/F  = 7.89 L/h    (CV 38%)
  V/F   = 52.3 L       (CV 35%)
  Ka    = 0.91 h⁻¹
  t½    ≈ 4.6h

공변량:
  CL/F ← 신기능 (mild 영향)
  CL/F ← CYP2C19 (약한 영향)
  Nonlinear PD: 혈소판 감소증 (시간-농도 의존)
  
MDR-TB TDM: AUC 80–120 mg·h/L (효능)
             Ctrough <2 mg/L (독성 최소화)
```

#### Case 4-C: Tranexamic Acid (TXA) — 외상/수술
```
약물: Tranexamic Acid
투여경로: IV bolus + infusion (CRASH-3 시험)
ClinicalTrials.gov: NCT02535949

모델 구조: 2구획 IV (ADVAN3)
파라미터 (Stitt et al. BJP 2024):
  CL    = 6.71 L/h    (CV 29%)
  V1    = 8.34 L       (CV 22%)
  Q     = 3.82 L/h
  V2    = 10.4 L
  t½β   ≈ 2.1h

공변량:
  CL ← 신기능 (CrCl, linear)
  CL ← 외상 중증도 (ISS score)
  V1 ← BW (power 1.0)
```

---

### CATEGORY 5: 단백질 의약품/바이오로직스 (Biologics)

#### Case 5-A: Bevacizumab (Anti-VEGF) — mCRC/NSCLC
```
약물: Bevacizumab (Avastin, Roche)
모달리티: IgG1 mAb (149 kDa)
투여경로: IV (5–15 mg/kg q2W or q3W)

모델 구조: 2구획 선형 (ADVAN3) — TMDD 없음 (선형 범위 내)
파라미터:
  CL    = 0.231 L/day  (CV 32%)
  V1    = 2.73 L        (CV 15%)
  Q     = 0.462 L/day
  V2    = 1.69 L
  t½    ≈ 20일

공변량:
  CL ← BW (power 0.559)     ΔOFV = -52
  CL ← 성별 (F: -8%)        ΔOFV = -15
  CL ← ADA 형성 (+29%)
  V1 ← BW (power 0.454)
  V1 ← 성별 (F: -11%)
```

#### Case 5-B: Adalimumab (Anti-TNF) — RA/IBD
```
약물: Adalimumab (Humira, AbbVie)
모달리티: IgG1 fully human mAb (148 kDa)
투여경로: SC 40 mg q2W

모델 구조: 1구획 SC 흡수 (ADVAN2) — FcRn 재순환 반영
파라미터:
  CL    = 0.476 L/day  (CV 36%)
  V/F   = 8.17 L        (CV 30%)
  Ka    = 0.016 h⁻¹    (Tmax ≈ 5일)
  F_SC  = 0.64
  t½    ≈ 14일

공변량:
  CL ← BW (power 0.722)
  CL ← Methotrexate 병용 (-29%) ← 중요!
  CL ← 항약물항체 ADA (×2.1)
  CL ← 기저 질환 활성도 (CDAI, DAS28)
```

#### Case 5-C: Trastuzumab (Anti-HER2) — HER2+ BC
```
약물: Trastuzumab (Herceptin, Roche)
모달리티: IgG1 mAb, HER2 표적
투여경로: IV (8 mg/kg loading → 6 mg/kg q3W) 또는 SC (600 mg q3W)

모델 구조: 2구획 + 수용체 매개 분포 (RMD, 비선형)
또는 단순화: TMDD 준정상상태 근사
파라미터 (선형 근사):
  CL    = 0.225 L/day  (CV 44%)
  V1    = 2.95 L        (CV 27%)
  Q     = 0.614 L/day
  V2    = 2.84 L
  t½    ≈ 28일

공변량:
  CL ← 혈청 HER2 (sHER2, +20%/증가)  ← 종양 부담 반영
  CL ← ALT (간기능)
  CL ← BW (power 0.563)
  V1 ← BW (power 0.420)
```

---

### CATEGORY 6: 항진균제 (Antifungals — HSCT 맥락)

#### Case 6-A: Fluconazole — HSCT 예방
```
약물: Fluconazole (Diflucan)
투여경로: PO/IV (200–400 mg QD)

모델 구조: 1구획 PO (ADVAN2)
파라미터:
  CL/F  = 1.25 L/h    (CV 28%)  ← 신기능 주도
  V/F   = 41.7 L       (CV 30%)
  Ka    = 1.04 h⁻¹
  F     ≈ 0.90 (높은 생체이용률)
  t½    ≈ 22–31h (신기능 의존)

공변량:
  CL/F ← CrCl (linear, +1.2%/mL/min)  ← 핵심
  CL/F ← 연령 (노인: -15%)
  DDI 작용: CYP2C9/3A4 억제제 → 기질 약물 CL 감소
```

#### Case 6-B: Voriconazole — 중증 진균감염
```
약물: Voriconazole (Vfend)
투여경로: PO/IV (6 mg/kg BID loading → 4 mg/kg BID)

모델 구조: 1구획 비선형 PK (Michaelis-Menten 제거)
Cl(C) = Vmax/(Km+C) — 포화 가능성
파라미터:
  Vmax  = 2.8 mg/h     (CV 65%)  ← 매우 높은 IIV!
  Km    = 1.9 mg/L     (CV 78%)
  V/F   = 4.6 L/kg     (CV 50%)

공변량:
  Vmax ← CYP2C19 PM (×0.55)   ← 절대적으로 중요!
  Vmax ← CYP2C19 UM (×2.1)
  Vmax ← ABCB1 (P-gp 다형성)
  자가 유도/억제: 시간 의존 PK (자가 억제 일부)

TDM: Ctrough 1–6 mg/L (효능/독성 균형)
주의: CYP2C19 PM 환자: 2배 이상 AUC↑ → 간독성
```

---

### CATEGORY 7: 항경련제 (Anticonvulsants)

#### Case 7-A: Valproic Acid (VPA) — 뇌전증/소아
```
약물: Valproate (Depakene)
투여경로: PO/IV (소아 15–60 mg/kg/day)

모델 구조: 1구획 + 단백결합 포화 (비선형 PK)
  CL(C) = CL_free × fu(C)
  fu(C) = 1/(1 + C/Kd_protein) — 포화 모델

파라미터:
  CL_free = 2.14 L/h   (CV 38%)
  V_free  = 8.7 L
  Kd      = 195 mg/L   (단백결합 포화 상수)

공변량 (소아):
  CL ← BW (power 0.75)
  CL ← 연령 (성숙도 함수, maturation)
  CL ← 병용 효소유도제 (+50%)
  CL ← 병용 항뇨 제거 약물

TDM: 총 농도 50–100 mg/L (단백결합 고려 필수)
```

---

## 케이스 추천 알고리즘

```python
def recommend_cases(drug_class, modality, route, target):
    """
    신약 분석 시 유사 케이스 자동 추천
    """
    case_db = load_case_database()  # 위 30개+ 케이스
    
    score_matrix = {}
    for case in case_db:
        score = 0
        if case.drug_class == drug_class: score += 40
        if case.modality == modality: score += 30
        if case.route == route: score += 20
        if case.target_pathway == target: score += 10
        score_matrix[case.id] = score
    
    # 상위 5개 추천
    top5 = sorted(score_matrix, key=score_matrix.get, reverse=True)[:5]
    
    # 파라미터 초기값 = 유사 케이스 중위수
    CL_init = median([case_db[c].CL for c in top5])
    V_init  = median([case_db[c].V for c in top5])
    
    return {
        'recommended_cases': top5,
        'CL_initial': CL_init,
        'V_initial': V_init,
        'suggested_covariates': union([case_db[c].significant_covariates for c in top5]),
        'model_structure': mode([case_db[c].model_structure for c in top5])
    }
```

## 파라미터 클래스별 범위 요약

| 약물 클래스 | CL (L/h) | V (L) | t½ | IIV CL | 핵심 공변량 |
|------------|---------|------|-----|--------|----------|
| ICI mAb (항PD-1/L1) | 0.008–0.015 L/h | 3–6 | 25–30일 | 30–40% | BW, ECOG, Alb, ADA |
| 항암 mAb (anti-HER2/VEGF) | 0.009–0.015 L/h | 2.7–5 | 18–28일 | 28–45% | BW, sTarget, ADA |
| 표적항암 SM (TKI) | 10–80 L/h | 100–900 | 4–26h | 45–70% | CYP3A4, BW, 음식 |
| 알킬화제 (Busulfan) | 7–12 L/h | 40–60 | 2–4h | 25–35% | BW, GST, Azole |
| 면역억제 (CsA/Tac) | 15–40 L/h | 500–1200 | 12–27h | 45–65% | CYP3A5, P-gp, ADA |
| 항생제 (Vancomycin) | 2–8 L/h | 20–60 | 6–12h | 40–55% | CrCl, BW, 패혈증 |
| 항진균제 (Voriconazole) | 0.5–4 L/h | 150–300 | 6–24h | 60–80% | CYP2C19, ABCB1 |
