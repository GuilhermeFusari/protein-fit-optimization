# Tarefa 1 — Efeito do número de inícios

**Pergunta:** com quantos inícios o ICP-SAXS deixa de ser significativamente
pior que o cifsup NSD em Chamfer?

**Resposta curta:** com **100** inícios, nas três sementes testadas. Com 3,
10 e 30 inícios ele continua significativamente pior (p < 0,05 em todas as
sementes). Com 100, o Chamfer médio (4,994–4,997 Å) é igual ao do cifsup NSD
(4,997 Å) e o Wilcoxon dá p = 0,066–0,088. Isso é "não significativamente
diferente", não "equivalente" (ver ressalvas).

Nada no pacote foi alterado; o default continua `--restarts 3`.

## Como foi medido

* `efeito_inicios.py`: roda `saxs_icp.align` com os parâmetros publicados
  nas 50 entradas, variando só `restarts` ∈ {3, 10, 30, 100}. Semente-base
  42 (a publicada) e, para medir a sensibilidade ao sorteio, 1042 e 2042.
* **Checagem:** 3 inícios + semente 42 reproduz `final_lam02.csv` com 0
  diferenças.
* cifsup NSD: valores lidos (sem modificar) de
  `benchmark/benchmark_cifsup_lam02.csv`.
* Wilcoxon pareado bilateral nas 50 entradas, como em
  `scripts/gerar_figuras_v2.py`. Com os CSVs publicados ele reproduz os
  p-valores do artigo (Chamfer 9,4e-6; fração fora 0,078; cobertura 0,001;
  separação 0,126).
* `analise_inicios.py` gera todas as tabelas abaixo (`analise_inicios.txt`).

## Resultado com a semente publicada (42)

| Inícios | Chamfer (Å) | p vs NSD | Fração fora | p | Cobertura | p | Sep. centroide (Å) | p |
|---|---|---|---|---|---|---|---|---|
| 3 (publicado) | 5,170 | 9,4e-6 | 0,124 | 0,078 | 0,528 | 0,001 | 4,733 | 0,126 |
| 10 | 5,053 | 3,2e-6 | 0,122 | 0,404 | 0,537 | 0,039 | 3,954 | 0,466 |
| 30 | 5,008 | 0,017 | 0,119 | 0,852 | 0,541 | 0,038 | 3,814 | 0,058 |
| 100 | 4,994 | 0,066 | 0,116 | 0,151 | 0,544 | 0,727 | 4,031 | 0,290 |
| cifsup NSD | 4,997 | | 0,123 | | 0,538 | | 4,197 | |

Também a cobertura deixa de ser significativamente pior com 100 inícios
(p = 0,73). Fração fora e separação já não eram significativas com 3.

## Robustez ao sorteio (Chamfer médio, p contra cifsup NSD)

| Inícios | semente 42 | semente 1042 | semente 2042 |
|---|---|---|---|
| 3 | 5,170 (p = 9,4e-6) | 5,203 (p = 2,4e-6) | 5,150 (p = 4,5e-6) |
| 10 | 5,053 (p = 3,2e-6) | 5,052 (p = 7,8e-4) | 5,025 (p = 1,7e-4) |
| 30 | 5,008 (p = 0,017) | 5,006 (p = 0,023) | 5,004 (p = 0,019) |
| 100 | 4,994 (p = 0,066) | 4,997 (p = 0,065) | 4,997 (p = 0,088) |

A conclusão é a mesma nas três sementes. Com 3 inícios o resultado depende
do sorteio (5,150–5,203 Å); com 100, quase nada (4,994–4,997 Å).

## O ganho vem da busca

* Mais inícios melhoram o Chamfer de forma significativa em relação a 3, em
  todas as sementes (100 inícios: −0,15 a −0,21 Å, p < 1e-7).
* O custo da Eq. 1 cai junto (média 7,13 → 6,90 com a semente 42): a busca
  encontra poses melhores pelo próprio objetivo, e o Chamfer acompanha.
  Com 3 inícios, parte das 50 poses publicadas não é o mínimo que o próprio
  método consegue achar.

## Ressalvas

1. **"Não significativo" não é "equivalente".** Com 100 inícios o cifsup NSD
   ainda tem Chamfer menor em 33 das 50 entradas (o ICP-SAXS em 17). A média
   empata porque o ICP-SAXS ganha por margens maiores nas entradas em que
   ganha. Para afirmar equivalência seria preciso um teste de equivalência
   (TOST) com uma margem definida antes, o que não fiz.
2. **Comparações múltiplas.** São 4 níveis × 3 sementes, sem correção. Com
   correção (Holm/Bonferroni), 30 inícios (p ≈ 0,02) provavelmente também
   deixaria de ser significativo. Escolher o nível de inícios depois de ver
   os p-valores é um ponto que um revisor pode levantar.
3. **Custo computacional.** 100 inícios custam ~33 vezes mais que 3 (≈ 6 min
   para as 50 entradas em 8 núcleos, contra ≈ 18 s).
4. **Inconsistência pequena na comparação publicada.** Em 9 das 50 entradas o
   envelope tem mais de 5000 pontos. Nelas o Chamfer do ICP-SAXS é medido
   contra o envelope amostrado a 5000 pontos, e o do cifsup contra o
   envelope completo (`benchmark_cifsup.py`). Medindo os dois contra o
   envelope completo, o Chamfer do ICP-SAXS muda ~−0,005 Å e os p-valores
   quase não mudam (3 inícios: 1,2e-5 em vez de 9,4e-6). Com 30 inícios,
   p = 0,033; com 100, p = 0,11. A conclusão não muda.

## Arquivos

* `efeito_inicios.py` — gera `efeito_inicios.csv` (600 linhas: 50 entradas
  × 4 níveis × 3 sementes). Leva ~25 min em 8 núcleos.
* `analise_inicios.py` — gera as tabelas (`analise_inicios.txt`).
* `efeito_inicios.log` — saída da execução, incluindo a checagem contra
  `final_lam02.csv`.
