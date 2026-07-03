---
name: pk-combo-therapy-agent
description: >
  약동학(PK) 모델링 기반 병용요법 설계 및 최적화 AI 에이전트 스킬.
  PopPK 프레임워크(1구획/2구획, FOCE, NONMEM 구조)를 활용하여 약물 간 PK 상호작용을 예측하고,
  시뮬레이션 기반으로 최적 병용 용량 및 일정을 제안함.

  다음 요청에 반드시 사용:
  - "병용요법 PK 분석", "약물 상호작용 예측", "DDI 시뮬레이션", "combination therapy PK"
  - "두 약물의 병용 용량 최적화", "병용 시 AUC 변화 예측"
  - "항암제 병용 PK", "면역억제제 병용", "항생제 병용 요법 설계"
  - "CYP 억제/유도 기반 DDI 예측", "GST 효소 병용 영향", "트랜스포터 매개 DDI"
  - "시너지 vs 拮抗 PK 평가", "병용 안전성 마진", "TDM 기반 병용 조정"
  - 약물 이름 두 개 이상 + "같이", "함께", "병용", "combination" 키워드
  - "PopPK 공변량으로 병용약물 추가", "공변량 모델에 DDI 항 추가"
  - "약동학 예측해서 병용요법 제안해줘" 류의 요청

  출력물:
  1. PK 파라미터 추정 테이블 (각 약물 단독 + 병용)
  2. DDI 분류 및 메커니즘 분석 보고서
  3. 시뮬레이션 기반 최적 용량 제안 (Interactive HTML 또는 xlsx)
  4. NONMEM 컨트롤 스트림 코드 스니펫
  5. 병용 안전성 평가 및 TDM 권고안

compatibility:
  - requires: bash_tool, create_file, present_files
  - optional: web_search (DDI 문헌 검색), PubMed (최신 PK 논문)
---

# PK Combo Therapy Agent

## 개요

이 스킬은 **Busulfan PopPK 분석 프레임워크**를 범용 약물 병용요법 PK 에이전트로 확장한 것입니다.
단일 약물 PopPK → 병용 DDI 모델링 → 시뮬레이션 → 최적 용량 제안까지 완전한 파이프라인을 제공합니다.

---

## 에이전트 실행 파이프라인 (9단계)

### Stage 0 — 약물 정보 파악
```
사용자 입력에서 추출:
  A. 약물 목록 (약물명 A, B, C...)
  B. 적응증 / 환자군 (종양학, 이식, 감염 등)
  C. 투여 경로 (IV / PO / SC / IM)
  D. 이미 알려진 PK 파라미터 (선행연구, 문헌)
  E. 우려되는 DDI 메커니즘 (CYP / GST / P-gp / UGT 등)

없으면: references/drug-pk-database.md의 기본값 적용
문헌 검색 필요 시: web_search + PubMed:search_articles 병렬 호출
```

### Stage 1 — 단독 약물 PK 프로파일 구축
```
각 약물에 대해:
  1. 구획 모델 결정 (1-CMT vs 2-CMT vs PBPK)
     → references/structural-model-selection.md 참조
  2. 핵심 PK 파라미터 추출:
     - CL (전신 청소율), V (분포용적)
     - ka (흡수속도, PO만), t1/2, Tmax, Cmax
     - 단백결합률, 주 대사효소, 배설경로
  3. IIV 정보 (CV% for CL, V)
  4. 특수 집단 조정 (신기능, 간기능, 연령, 체중)
```

### Stage 2 — DDI 메커니즘 분류
```
DDI 유형 분류 (references/ddi-classification.md):

[약물대사 DDI]
  CYP 억제: 가역적(경쟁적) / 시간의존적(TDI) / 비가역적(MI)
  CYP 유도: PXR/CAR 활성화 → mRNA/단백 증가 (수일 지연)
  UGT 억제/유도
  FMO 기질 충돌

[트랜스포터 DDI]
  P-gp (MDR1): 뇌, 장관, 신세뇨관
  BCRP: 장관 흡수, BBB
  OATP1B1/1B3: 간 uptake
  OCT2/MATE: 신 배설

[약력학적 상호작용]
  QT 연장 중첩
  혈액독성 중첩 (골수억제)
  신독성, 간독성 중첩

[단백결합 경쟁]
  Albumin 결합 경쟁 → 유리형 증가

분류 후 → DDI 심각도 등급 결정:
  Grade 1 (Minor): ΔAUC <25%
  Grade 2 (Moderate): ΔAUC 25~100%
  Grade 3 (Major): ΔAUC >100% OR 임상적으로 위험
  Grade 4 (Contraindicated): 절대 병용 금기
```

### Stage 3 — PopPK DDI 모델 구조 설계
```
DDI를 공변량으로 모델링:

[CL에 대한 DDI 공변량 항]
  억제제 병용: CL_A = TVCL_A × exp(θ_inh × I_DrugB)
               θ_inh < 0 (CL 감소)
  유도제 병용: CL_A = TVCL_A × exp(θ_ind × I_DrugB)
               θ_ind > 0 (CL 증가)
  복합 DDI:   CL_A = TVCL_A × exp(θ_inh × I_B) × exp(θ_ind × I_C)

[Michaelis-Menten 억제 모델 (in vitro → in vivo scaling)]
  R = 1 + (I_u / Ki)              ← 가역적 억제
  R_TDI = 1 + (kinact/kdeg)×(I_u/(I_u+KI))  ← TDI

[흡수 DDI (PO 약물)]
  ka_ratio = ka_alone × (1 - Emax_DDI × C_B / (EC50 + C_B))

NONMEM 코드 스니펫 → references/nonmem-ddi-templates.md
```

### Stage 4 — PK 시뮬레이션
```
시뮬레이션 설정:
  N = 1,000 가상 환자 (IIV 적용)
  공변량 분포: 목표 환자군 특성 반영
  용량 범위: 0.5× ~ 2× 표준 용량 스윕

산출 지표:
  - AUC0-∞, AUC0-τ, Cmax, Tmax, Ctrough
  - DDI ratio: AUC_combo / AUC_alone (목표: 0.5~2.0)
  - Therapeutic Index 여유: (독성 AUC - 예측 AUC) / 예측 AUC
  - 목표 노출 달성률 (% patients within target range)

시뮬레이션 코드 → references/simulation-engine.md (Python/R)
```

### Stage 5 — 시너지/길항 평가 (PD 통합)
```
Loewe Additivity 모델:
  d_A/D_A + d_B/D_B = 1 (additivity)
  < 1: 시너지, > 1: 길항

Bliss Independence 모델:
  E_combo = E_A + E_B - E_A × E_B

CI (Combination Index) 계산:
  CI < 1: 시너지
  CI = 1: 가산
  CI > 1: 길항

→ PK/PD 통합 시뮬레이션으로 최적 병용 비율 탐색
```

### Stage 6 — 안전성 마진 및 TDM 권고
```
각 약물별:
  Safety margin = (NOAEL 또는 독성 threshold) / 예측 Cmax_combo
  목표: Safety margin ≥ 2.0

TDM 권고:
  채혈 시점: 가장 중요한 PK 파라미터를 포착하는 sparse 설계
  모니터링 주기: DDI 심각도에 따라 조정
  용량 조정 알고리즘: AUC-based 또는 Ctrough-based

특수 집단 권고:
  신기능 저하 (eGFR <30): 용량 및 간격 조정
  간기능 저하 (Child-Pugh B/C): CYP 기질 약물 감량
  노인 (≥65세): 낮은 초기 용량, 느린 증량
```

### Stage 7 — 최적 병용 프로토콜 제안
```
출력 형식 결정 (사용자 요청에 따라):
  Option A: Interactive HTML 대시보드 (시뮬레이션 슬라이더 포함)
  Option B: xlsx 보고서 (파라미터 테이블 + VPC 플롯)
  Option C: PDF 요약 보고서
  Option D: NONMEM 코드 + 분석 스크립트

최적 병용 프로토콜 포함 내용:
  1. 투여 순서 (sequence): 어떤 약물을 먼저?
  2. 투여 간격 (interval): 동시투여 vs 시차 투여
  3. 용량 조정 (dose adjustment): DDI 보정 용량
  4. 모니터링 계획 (monitoring plan): TDM 일정
  5. 중단 기준 (stopping criteria): 독성 임계값
```

### Stage 8 — 출력물 생성
```
참조: references/output-templates.md

[Interactive HTML 대시보드]
  - Anthropic API 활용한 AI-powered PK 시뮬레이터
  - 슬라이더: 용량, 체중, 신기능, DDI 억제 강도
  - 실시간 농도-시간 곡선 업데이트
  - AUC / Cmax / Ctrough 자동 계산
  - 목표범위 달성 여부 시각 표시

[NONMEM 컨트롤 스트림]
  - DDI 공변량 모델 완성본
  - $PK, $ERROR, $THETA, $OMEGA, $SIGMA 블록

[임상 권고 보고서]
  - Executive summary (1페이지)
  - DDI 메커니즘 + 예측 크기
  - 용량 조정 표
  - TDM 프로토콜
  - 참고문헌 (PubMed + FDA 가이드라인)
```

### Stage 9 — 검증 및 불확실성 보고
```
모델 불확실성 정량화:
  - In vitro → in vivo 스케일링 오차 (통상 2~5배)
  - IIV가 DDI 크기에 미치는 영향
  - 환자군 특성의 불확실성

검증 권고:
  - 외부 데이터셋으로 VPC 수행
  - Bootstrap CI 보고
  - 민감도 분석 (Ki 값 ±50%)

한계 명시:
  - 생리학적 복잡성 단순화
  - 임상 검증 필요성 명시
  - 규제 제출용으로는 추가 연구 필요
```

---

## 약물 DB 조회 우선순위

1. **사용자 제공 데이터** (임상/전임상 결과) — 최우선
2. **references/drug-pk-database.md** — 주요 약물 PK 파라미터
3. **PubMed:search_articles** — 최신 문헌
4. **web_search** — FDA label, 제품 정보
5. **추정값 사용 시** → 반드시 불확실성 명시

---

## 지원 약물 클래스 및 주요 DDI 패턴

| 클래스 | 대표 약물 | 주요 DDI 메커니즘 |
|--------|---------|--------------|
| 알킬화제 | Busulfan, Cyclophosphamide | GST (GSTA1), CYP3A4 |
| 칼시뉴린 억제제 | Cyclosporine, Tacrolimus | CYP3A4 S/I, P-gp |
| 항진균제 | Fluconazole, Voriconazole | CYP2C9/3A4 강력 억제 |
| 항경련제 | Phenytoin, Carbamazepine | CYP3A4 유도 |
| 항생제 | Rifampicin, Linezolid | CYP3A4 강력 유도 / MAO |
| 표적항암제 | Imatinib, Ibrutinib | CYP3A4, P-gp |
| mTOR 억제제 | Sirolimus, Everolimus | CYP3A4, P-gp 기질 |
| 항바이러스제 | Ritonavir, Cobicistat | CYP3A4 억제 (부스터) |

---

## 에이전트 응답 형식 규칙

1. **항상 DDI 심각도 등급 먼저 표시** (Grade 1~4)
2. **수치 예측에는 신뢰구간 포함** (예: AUC ratio 1.8 [1.3–2.5])
3. **임상적 행동 권고를 명확히** (용량 X% 감량, 모니터링 강화 등)
4. **면책 조항 포함**: "이 예측은 시뮬레이션 기반이며 임상 검증 필요"
5. **참고 문헌 인용**: FDA DDI 가이드라인 (2020), EMA DDI 가이드라인 (2012)

---

## 핵심 참조 파일

```
references/
├── drug-pk-database.md          ← 주요 약물 PK 파라미터 DB
├── ddi-classification.md        ← DDI 유형 분류 체계
├── structural-model-selection.md ← 구획 모델 선택 기준
├── nonmem-ddi-templates.md      ← NONMEM DDI 공변량 코드
├── simulation-engine.md         ← Python PK 시뮬레이션 엔진
├── output-templates.md          ← HTML/xlsx 출력 템플릿
└── fda-ddi-guidance-2020.md     ← FDA In Vitro DDI Studies Guidance
```
