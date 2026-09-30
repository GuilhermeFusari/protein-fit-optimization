# Tarefa 2 — Pesos do Procrustes vs Eq. 1

**Resposta curta:** o parecer está certo em que os pesos do Procrustes não
seguem a Eq. 1. O **custo** avaliado a cada iteração é exatamente a Eq. 1,
mas os **pesos** usados para calcular a rotação e a translação não. Uma
versão com os pesos da Eq. 1 (`eq1_ls`) otimiza melhor a própria Eq. 1
(custo menor em 42/50 entradas com 30 inícios, p = 3,8e-7), mas **não muda o
Chamfer de forma significativa** (3 inícios: 5,202 vs 5,170 Å, p = 0,28;
30 inícios: 4,987 vs 5,008 Å, p = 0,17). A versão IRLS, que minimiza a Eq. 1
como média de distâncias, saiu pior.

Nada no pacote foi alterado. A decisão de trocar é sua.

## Premissa: qual é a Eq. 1

O artigo (.docx) não está nesta máquina. Usei a equação do README (seção
"Choosing λ"), que é a mesma implementada em `author_cost`:

    J = mean_e d(e, P) + λ · mean_p d(p, E)       λ = 0,2
        (preenchimento)    (vazamento)

Se a Eq. 1 do artigo for diferente, esta análise precisa ser refeita.

## Onde o código diverge da Eq. 1

Arquivo `src/saxs_icp/core.py`, função `icp`, modo `author` (idêntico ao
`scripts/icp_saxs_v2.py`):

| Linha | O que faz | Coerente com a Eq. 1? |
|---|---|---|
| 117–118 | `cost = fill + penalty * leak` para escolher a melhor iteração | **Sim**, é a Eq. 1. |
| 129 | peso de cada par de vazamento (átomo → envelope): `w_i = 1 + λ · d_i / mean(d)` | **Não.** |
| 140–145 | pares de preenchimento (envelope → átomo) com peso `w_j = d_j / mean(d_env)` | **Não.** |
| 161 | `wA = wA / wA.sum()`: normaliza os dois blocos juntos | **Não.** |

Três divergências concretas:

1. **Peso relativo dos dois termos.** Na Eq. 1, vazamento/preenchimento =
   λ = 0,2. No código, a soma dos pesos do bloco de vazamento é
   N_P(1 + λ) e a do preenchimento é N_E. A razão efetiva é
   N_P(1 + λ)/N_E e depende do número de pontos: no benchmark, mediana 0,65,
   de 0,06 a 3,46 entre entradas.
2. **O "1 +" no vazamento.** Mesmo com λ pequeno, o vazamento tem peso ≥ 1
   por par. Na Eq. 1 ele pesa λ.
3. **Peso crescente com a distância.** Dentro de cada bloco, o peso é
   proporcional à distância (d/mean). Minimizar Σ w·d² com w ∝ d aproxima
   Σ d³, que dá ainda mais importância aos pontos distantes. A Eq. 1 é uma
   média de distâncias (Σ d), cuja reponderação coerente é w ∝ 1/d, o
   oposto.

Consequência prática: o passo de Procrustes minimiza outra função, e o
custo da Eq. 1 só entra na escolha da melhor iteração. O docstring do `icp`
("os pesos entram no passo de Procrustes") está certo em que existem pesos,
mas eles não são os da Eq. 1.

## As versões testadas (`pesos_procrustes.py`)

Tudo idêntico ao método publicado (custo J para escolher a iteração, 50
iterações, inícios, enantiômeros, amostragem, sementes); muda só o peso de
cada par no Procrustes:

| Variante | Vazamento (por par) | Preenchimento (por par) |
|---|---|---|
| `atual` | 1 + λ·d/mean(d) | d/mean(d_env) |
| `eq1_ls` | λ / N_P | 1 / N_E |
| `eq1_irls` | λ / (N_P · max(d, 0,1 Å)) | 1 / (N_E · max(d, 0,1 Å)) |

`eq1_ls` aplica os pesos da Eq. 1 à forma quadrática usual do ICP.
`eq1_irls` é a reponderação para a Eq. 1 como média de distâncias (tipo
Weiszfeld); o ε = 0,1 Å evita pesos infinitos.

**Checagem:** `atual` com 3 inícios reproduz `final_lam02.csv` com 0
diferenças.

## Resultados (50 entradas, Wilcoxon pareado)

**3 inícios (configuração publicada)**

| Variante | Chamfer (Å) | Custo J | Fração fora | Cobertura | Sep. centroide (Å) | Chamfer: p vs atual | J: p vs atual | J menor que atual | Chamfer: p vs cifsup NSD |
|---|---|---|---|---|---|---|---|---|---|
| atual | 5,170 | 7,131 | 0,124 | 0,528 | 4,733 | — | — | — | 9,4e-6 |
| eq1_ls | 5,202 | 7,125 | 0,127 | 0,531 | 4,838 | 0,28 | 0,023 | 35/50 | 1,1e-8 |
| eq1_irls | 5,406 | 7,396 | 0,142 | 0,533 | 6,550 | 1,7e-5 | 6,6e-5 | 15/50 | 3,9e-10 |

**30 inícios**

| Variante | Chamfer (Å) | Custo J | Fração fora | Cobertura | Sep. centroide (Å) | Chamfer: p vs atual | J: p vs atual | J menor que atual | Chamfer: p vs cifsup NSD |
|---|---|---|---|---|---|---|---|---|---|
| atual | 5,008 | 6,920 | 0,119 | 0,541 | 3,814 | — | — | — | 0,017 |
| eq1_ls | 4,987 | 6,855 | 0,118 | 0,546 | 4,284 | 0,17 | 3,8e-7 | 42/50 | 0,105 |
| eq1_irls | 5,087 | 6,984 | 0,128 | 0,551 | 5,169 | 0,030 | 0,19 | 21/50 | 0,002 |

cifsup NSD: Chamfer médio 4,997 Å.

## Interpretação

1. **`eq1_ls` otimiza melhor o objetivo declarado.** Com pesos coerentes, o
   custo da Eq. 1 fica menor na maioria das entradas (35/50 com 3 inícios,
   42/50 com 30). Isso é o esperado quando o passo de Procrustes passa a
   descer a mesma função que é avaliada.
2. **Mas o Chamfer não muda significativamente.** A diferença de Chamfer
   entre `eq1_ls` e `atual` não é significativa com 3 nem com 30 inícios. A
   separação de centroides piora um pouco com `eq1_ls` (30 inícios: 4,28 vs
   3,81 Å; não testei a significância).
3. **Com 30 inícios, `eq1_ls` deixa de ser significativamente pior que o
   cifsup NSD** (p = 0,105; a versão atual dá p = 0,017). É um indício, não
   uma conclusão: são comparações múltiplas e a diferença para `atual` não é
   significativa.
4. **`eq1_irls` é pior** em Chamfer e em J com 3 inícios. Não investiguei a
   causa (ε, interação com a troca de correspondências a cada iteração).
   Como está, não é candidata.

**Leitura para o parecer:** a divergência existe e agora está descrita
exatamente. Corrigir os pesos torna o método coerente com a equação do
artigo sem piorar o resultado. Com 3 inícios, porém, a busca pesa mais que
os pesos (Tarefa 1). Se o método for trocado, o caminho mais defensável
parece ser `eq1_ls` junto com mais inícios, mas isso exige refazer o
benchmark, a calibração de λ (o ótimo de λ pode mudar com pesos
diferentes) e as figuras. A decisão é sua.

## Arquivos

* `pesos_procrustes.py` — as três variantes. Leva ~8 min em 8 núcleos.
* `pesos_procrustes.csv` — 300 linhas (50 entradas × 3 variantes × 2
  níveis de inícios).
* `analise_procrustes.py` → `analise_procrustes.txt` — as tabelas acima.
* `pesos_procrustes.log` — saída, com a checagem contra `final_lam02.csv`.
