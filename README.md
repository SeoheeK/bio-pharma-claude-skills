# Bio-Pharma Claude Skills

<img src="assets/banner.svg" alt="Bio-Pharma Claude Skills" width="100%" />

**AI 기반 신약개발 · 임상약동학(PK) · 규제과학 업무를 위한 Claude Skill 패키지 모음**

타겟 발굴부터 IND 신청, 임상시험 설계, 이득-위험 심사보고서 작성까지 — 바이오/제약 R&D 파이프라인의 각 단계를 담당하는 8개의 독립 Claude Skill을 제공합니다. 각 스킬은 [Claude Skills](https://docs.claude.com) 표준 형식(`SKILL.md` + 참조 문서)으로 작성되어, Claude.ai / Claude Code / API에 그대로 업로드해 사용할 수 있습니다.

> 이 저장소는 KRIBB 국가생명공학정책연구센터에서 바이오 정책·AI 신약개발 리서치를 수행하며 축적한 도메인 지식을 스킬 형태로 구조화한 것입니다.

---

## 목차

| # | 스킬 | 한 줄 요약 |
|---|---|---|
| 01 | [`ai-drug-discovery-pipeline`](#01--ai-drug-discovery-pipeline) | 타겟 선정 → AI 분자생성 → 11-Layer 스크리닝 → TPP → DBTL → 의사결정 |
| 02 | [`benefit-risk-assessment`](#02--benefit-risk-assessment) | ICH M4E(R2)·CTD 2.5.6 기반 이득-위험 심사보고서 자동 작성 |
| 03 | [`biopharma-gng-decision`](#03--biopharma-gng-decision) | 개발 단계별 Go / Conditional Go / No-Go 스코어카드 |
| 04 | [`pk-allometry-ensemble`](#04--pk-allometry-ensemble) | 동물 PK → 인간 CL/Vd 예측 (10종 Allometry 앙상블) |
| 05 | [`pk-biomarker-db`](#05--pk-biomarker-db) | 타겟/MoA 기반 PD 바이오마커 추천 + PK/PD 모델 구조 |
| 06 | [`pk-case-study-db`](#06--pk-case-study-db) | 30+ 임상시험 PopPK 케이스 매칭 + 파라미터 초기값 제안 |
| 07 | [`pk-combo-therapy-agent`](#07--pk-combo-therapy-agent) | 병용요법 DDI 예측 + 최적 용량·일정 시뮬레이션 |
| 08 | [`pk-study-design-simulator`](#08--pk-study-design-simulator) | Monte Carlo 기반 PK 임상시험 설계 최적화 |

---

## 01 · `ai-drug-discovery-pipeline`

<img src="assets/ai-drug-discovery-pipeline.svg" alt="AI Drug Discovery Pipeline" width="100%" />

AI 기반 신약 발굴을 위한 **5단계 통합 파이프라인 프레임워크**입니다. 타겟·질환 유효성 확립(TVS)부터 Evo2·AlphaFold3·REINVENT를 활용한 분자 생성, 11-Layer 깔때기형 Virtual Screening, TPP(Target Product Profile) 수립, Bocheming+SDL 기반 DBTL 자동화 검증 루프, 그리고 임상·시장성·ESG를 종합한 Go/No-Go 의사결정까지 전 과정을 다룹니다.

- **트리거 예시**: "신약 파이프라인 설계해줘", "가상 스크리닝 필터 기준 알려줘", "REINVENT 보상 함수 어떻게 짜?", "TPP 수립해줘", "타겟 유효성 검증"
- **구성**: 5단계 개요 + 8개 상세 참조 문서(references/00~07) — 임상 바이오마커 군집분석, 타겟·질환 유효성, AI 생성 도구, 스크리닝 필터, DBTL 설계, 번역·의사결정, TPP 프레임워크, 필터 근거 논문 DB
- **연동 스킬**: `biopharma-gng-decision`(Go/No-Go), `pk-allometry-ensemble`(1상 시작용량), `pk-biomarker-db`(PD 마커), `bio-regulatory-analysis`(규제 경로)

---

## 02 · `benefit-risk-assessment`

<img src="assets/benefit-risk-assessment.svg" alt="Benefit-Risk Assessment" width="100%" />

**ICH M4E(R2)·식약처 CTD 2.5.6 공식 서식**과 **FDA 구조화 이득-위험 프레임워크(BRF)**에 기반해, 임상·전임상 유효성·안전성 데이터를 정량(NNT/NNH, MCDA, QALY, 안전역)·정성적으로 종합 평가하고 심사보고서를 작성합니다. 전임상 단계는 세포·오가노이드·마우스·랫드·토끼·개·원숭이 등 시험계별로 데이터를 구분해 평가합니다.

- **트리거 예시**: "이득-위험 평가해줘", "NNT NNH 계산", "MCDA 분석", "이 약물 허가해도 되나", "전임상 안전역 계산", "REMS 제안"
- **출력물**: 이득-위험 종합 심사보고서(식약처 CTD 2.5.6 구조, Markdown/.docx), FDA Figure 1 매트릭스, NNT/NNH/MCDA/QALY 정량분석, 참고문헌 부록
- **구분**: 개발단계 Go/No-Go 판단은 `biopharma-gng-decision`, 이 스킬은 "허가 여부·이득-위험 상회 여부"를 심사관 관점에서 평가

---

## 03 · `biopharma-gng-decision`

<img src="assets/biopharma-gng-decision.svg" alt="Biopharma Go/No-Go Decision" width="100%" />

전임상 → IND → Ph1 → Ph2 → Ph3 → NDA/BLA, 개발 단계별 기준 변수(Efficacy·Safety·CMC·Regulatory·rNPV)를 점수화하여 **Go / Conditional Go / No-Go**를 자동 산출하는 의사결정 알고리즘 스킬입니다.

- **트리거 예시**: "Go/No-Go 의사결정", "IND 신청 가능 여부 판단", "Ph2→Ph3 진입 판단", "Kill criteria 정리해줘", "이 데이터로 계속 개발해야 하나?"
- **출력물**: Go/No-Go 스코어카드 엑셀(.xlsx, 단계별 탭 + 판정 결과), Mermaid 의사결정 트리 플로우차트

---

## 04 · `pk-allometry-ensemble`

<img src="assets/pk-allometry-ensemble.svg" alt="PK Allometry Ensemble" width="100%" />

동물(Mouse·Rat·Rabbit·Dog·Monkey) → 인간 약동학 파라미터 스케일링을 위한 **Allometry 앙상블 스킬**. ChEMBL API에서 종간 PK 데이터를 자동 수집하고, Simple Allometry·Rule of Exponents·Brain-Weight 보정·MLP 보정 등 10가지 방법을 앙상블하여 인간 CL/Vd를 예측합니다. UniPlot(종간 Log-Log 시각화)을 함께 제공합니다.

- **트리거 예시**: "동물 PK → 인간 예측", "rule of exponents", "전임상 → FIH 스케일링", "ChEMBL PK 데이터 가져와줘", "rat dog monkey → human CL 예측"
- **출력물**: ChEMBL 종간 PK 데이터(parquet/CSV), 10종 방법별 예측값 + 앙상블 결과, UniPlot 시각화, Snakemake 워크플로우 자동 실행

---

## 05 · `pk-biomarker-db`

<img src="assets/pk-biomarker-db.svg" alt="PK Biomarker DB" width="100%" />

약물 표적(target)·작용기전(MoA)·치료 영역을 입력하면 활성/저해가 예측되는 **약동력학적(PD) 바이오마커를 추천**하고 PK/PD 모델 구조(직접·간접 반응 모델)를 제안합니다. ClinicalTrials.gov 프로토콜 + FDA 심사문서 + 주요 논문을 근거로 합니다.

- **트리거 예시**: "이 약물의 PD 마커 추천해줘", "target engagement 마커", "면역항암제 바이오마커", "독성 바이오마커 추천", "phase 1 바이오마커 전략"
- **출력물**: 약물/표적별 권장 PD 바이오마커(활성/저해 분류), PK/PD 모델 구조 제안, 마커별 측정법 + 채혈 시점, 단계별 바이오마커 전략, NONMEM 코드 스니펫

---

## 06 · `pk-case-study-db`

<img src="assets/pk-case-study-db.svg" alt="PK Case Study DB" width="100%" />

ClinicalTrials.gov 및 FDA 심사 데이터 기반, **30개 이상 실제 임상시험 약물의 PopPK 케이스 스터디 DB**. 새 약물 분석 시 유사 케이스를 자동 추천하고 파라미터 초기값(THETA/OMEGA/SIGMA)과 NONMEM 코드를 생성합니다.

- **트리거 예시**: "이 약물과 비슷한 PK 모델 있어?", "파라미터 초기값 추천", "mAb PK 모델", "항암제 PopPK 사례", "이 약물 클래스 대표 파라미터 알려줘"
- **출력물**: 유사 케이스 3~5개 자동 추천(약물명·NCT번호·모델 구조), 추천 초기 파라미터 세트, 문헌 기반 공변량 선택 근거, NONMEM 컨트롤 스트림 템플릿

---

## 07 · `pk-combo-therapy-agent`

<img src="assets/pk-combo-therapy-agent.svg" alt="PK Combo Therapy Agent" width="100%" />

**약동학(PK) 모델링 기반 병용요법 설계·최적화 에이전트.** PopPK 프레임워크(1/2-구획, FOCE, NONMEM)로 약물 간 CYP 억제/유도·트랜스포터 매개 DDI를 예측하고, 시뮬레이션으로 최적 병용 용량·일정을 제안합니다.

- **트리거 예시**: "병용요법 PK 분석", "DDI 시뮬레이션", "두 약물의 병용 용량 최적화", "항암제 병용 PK", "TDM 기반 병용 조정"
- **출력물**: 단독+병용 PK 파라미터 추정 테이블, DDI 분류·메커니즘 분석 보고서, 최적 용량 제안(Interactive HTML/xlsx), NONMEM 코드 스니펫, 병용 안전성·TDM 권고안

---

## 08 · `pk-study-design-simulator`

<img src="assets/pk-study-design-simulator.svg" alt="PK Study Design Simulator" width="100%" />

약물 특성(PK 파라미터, 전임상 데이터)을 입력하면 **최적 채혈설계(sparse/dense), 필요 환자 수, 용량 범위, 투여 일정**을 Monte Carlo 시뮬레이션으로 탐색해 최적 설계안 3가지를 비교 제안하는 PK 임상시험 설계 시뮬레이터입니다. Interactive HTML 대시보드로 예측 정밀도(RSE%)·파워·목표 AUC 달성률을 실시간 확인할 수 있습니다.

- **트리거 예시**: "PK 연구 설계", "몇 명 환자가 필요한가", "sparse sampling 설계", "First-in-human PK 설계", "FIH dose escalation 설계"
- **출력물**: Interactive HTML 연구설계 대시보드, 최적 설계안 3개 비교표, Power calculation 요약, FDA/ICH E14·EMA 가이드라인 준수 체크리스트

---

## 사용 방법

1. 원하는 스킬의 `.zip`(또는 `.md`) 파일을 다운로드합니다.
2. Claude.ai → **설정 → Capabilities → Skills**(또는 Claude Code `/skills`)에서 업로드합니다.
3. 스킬 설명(description)에 명시된 트리거 키워드로 대화를 시작하면 자동으로 활성화됩니다.

여러 스킬은 서로 연동되도록 설계되어 있습니다 (예: `ai-drug-discovery-pipeline` → `biopharma-gng-decision` → `benefit-risk-assessment` 순으로 Discovery부터 허가 심사까지 이어짐).

## 폴더 구조

```
bio-pharma-claude-skills/
├── README.md
├── assets/                          # README 다이어그램 (SVG)
├── ai-drug-discovery-pipeline.zip
├── benefit-risk-assessment.zip
├── biopharma-gng-decision.zip
├── pk-allometry-ensemble.md
├── pk-biomarker-db.md
├── pk-case-study-db.md
├── pk-combo-therapy-agent.md
└── pk-study-design-simulator.md
```

## License

이 저장소의 스킬은 리서치·참조 목적으로 공개되며, 실제 임상·규제 의사결정에는 반드시 전문가 검토를 거쳐야 합니다.
