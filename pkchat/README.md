# pkchat — 챗봇 기반 약동학(PK) 모델링 엔진

`bio-pharma-claude-skills`의 **문서형 PK 스킬**(`pk-*.md`)에 대응하는 **실제 실행 코드 레이어**입니다.
스킬이 "무엇을·왜"의 도메인 지식을 제공한다면, `pkchat`은 그것을 **통계/수치/ML 코드**로 구현하고
**챗봇(자연어) 인터페이스**로 구동합니다.

> 설계 원칙: **수치 코어는 Python 표준 라이브러리만으로 동작**합니다(numpy/scipy 불필요).
> 재현성·이식성이 높고, 어떤 환경에서도 그대로 실행됩니다. Claude API(`anthropic`)는
> 대화 레이어에만 **선택적으로** 사용됩니다.

---

## 아키텍처

```
자연어 요청
    │
    ▼
┌──────────────────────────────────────────────┐
│  chat/  — 챗봇 레이어                          │
│   • Claude 백엔드 (anthropic tool-use 루프)    │  ← ANTHROPIC_API_KEY 있으면 자동
│   • Local 백엔드 (규칙 기반 인텐트 라우터)      │  ← 의존성 0, 항상 동작
│   • tools.py — 6개 도구 (Anthropic 스키마 호환)│
└───────────────┬──────────────────────────────┘
                ▼   동일한 도구 표면(tool surface)
┌──────────────────────────────────────────────┐
│  수치 엔진 (pure Python)                        │
│   models / regimen / simulate  구획모델·ODE·시뮬 │
│   nca            비구획 분석 (AUC, t½, CL/Vz)    │
│   estimation     개별·2단계 집단 PK 파라미터 추정 │
│   optim          Nelder–Mead 최적화             │
│   ml/allometry   동물→인간 스케일링 앙상블(+ML)  │
│   casedb         PopPK 케이스 DB + 유사도 추천   │
└──────────────────────────────────────────────┘
```

두 챗봇 백엔드는 **같은 도구 레지스트리**(`chat/tools.py`)를 호출하므로 수치 결과는 백엔드와
무관하게 동일합니다. Claude가 있으면 언어 이해가 유연해지고, 없으면 로컬 라우터가 문서화된
요청 패턴을 처리합니다.

---

## 빠른 시작

```bash
# 대화형 REPL (백엔드 자동 선택)
python -m pkchat.cli

# 원샷 실행
python -m pkchat.cli "simulate CL=5 V1=30 dose=100 iv"

# 로컬 라우터 강제 (Claude 미사용)
python -m pkchat.cli --no-claude "scale to human: rat 0.25 3.5, dog 10 45, monkey 5 28"

# 전체 데모
python examples/demo.py

# 테스트 (22개, 표준 라이브러리만)
python -m unittest discover -s tests
```

Claude 백엔드를 쓰려면: `pip install anthropic` 후 `ANTHROPIC_API_KEY` 설정.
기본 모델은 `claude-opus-4-8`이며 `--model`로 변경 가능합니다.

---

## 구현된 기능

| 모듈 | 기능 | 핵심 |
|------|------|------|
| `models` | 1/2-구획 모델 (IV bolus/infusion/oral), 선형·비선형(Michaelis–Menten) | ODE(RK4) + 해석해(closed-form) 이중 구현, 상호 검증 |
| `simulate` | 임의 투약 요법(다회·병용 투여) 시뮬레이션, 집단 시뮬레이션 | 로그정규 IIV + 비례 잔차오차, 예측구간(percentile band) |
| `nca` | 비구획 분석 | linear-up/log-down AUC, 말단 반감기(adjusted-R² 창 선택), CL/Vz |
| `estimation` | 파라미터 추정 | 로그공간 가중최소제곱, 근사 RSE(%), 2단계 NLME(THETA/OMEGA 요약) |
| `ml/allometry` | 동물→인간 스케일링 | Simple/Fixed-b/Rule of Exponents/MLP/뇌중량/2종/bagged 회귀 앙상블, LOO |
| `casedb` | PopPK 케이스 DB | 15개 실약물, 클래스·모달리티·경로 유사도 추천 + 초기 파라미터 제안 |

### 정확도 검증
ODE 적분기는 1·2-구획 및 경구 모델에서 해석해와 **오차 ~1e-11**로 일치하며
(`tests/test_pkchat.py`), 추정기는 무잡음 데이터에서 참값 CL/V를 **오차 <0.1%**로 복원합니다.

---

## Python API 예시

```python
from pkchat import PKParameters, DosingRegimen, simulate, nca, fit_individual

# 2-구획 경구, 다회 투여 시뮬레이션
p = PKParameters(CL=4.0, V1=8.0, Q=3.0, V2=10.0, ka=1.2, F=0.9)
reg = DosingRegimen.multiple(200, interval=12, n_doses=6, route="oral")
res = simulate(p, reg, t_end=96)
print(res.cmax(), res.tmax())

# 관측 데이터로 CL/V 추정
fit = fit_individual(times, concentrations, DosingRegimen.single(100, "iv_bolus"),
                     estimate=("CL", "V1"))
print(fit.summary())
```

```python
from pkchat.chat import PKChatAgent
agent = PKChatAgent()            # Claude 있으면 사용, 없으면 로컬
print(agent.chat("simulate a 100 mg IV bolus with CL=5 and V=30"))
```

---

## 챗봇 도구 표면 (6개)

`simulate_pk`, `run_nca`, `fit_pk_model`, `allometric_scaling`,
`recommend_pk_case`, `simulate_population` — 모두 Anthropic tool-use JSON 스키마로
정의되어 있어(`pkchat.chat.tool_schemas()`) 그대로 Claude Messages API의 `tools`로 전달됩니다.

---

## 기존 스킬과의 관계

| 스킬(문서) | 대응 코드 |
|-----------|----------|
| `pk-allometry-ensemble` | `ml/allometry.py` (앙상블 실구현) |
| `pk-case-study-db` | `casedb.py` (구조화 DB + 추천 알고리즘) |
| `pk-combo-therapy-agent` | `simulate.py` 다회·병용 투여 + `estimation.py` |
| `pk-study-design-simulator` | `simulate_population` (Monte-Carlo, 예측구간) |

스킬은 도메인 지식과 NONMEM 템플릿을, `pkchat`은 규제 제출 전 **탐색·프로토타이핑용
실행 엔진**을 제공합니다.

---

## 한계 및 주의

- 2단계 NLME는 투명한 개별 피팅 기반 요약이며, **완전한 FOCE 혼합효과 추정 엔진이 아닙니다**.
  규제 제출용 최종 분석에는 NONMEM/Monolix 등을 사용하십시오.
- 모든 예측은 **모델 기반**이며 임상·규제 의사결정 전 반드시 전문가 검토가 필요합니다.
- 단위는 호출자가 일관되게 유지해야 합니다(관례: CL [L/h], V [L], dose [mg], time [h]).
