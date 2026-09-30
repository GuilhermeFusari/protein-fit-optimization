# Correções do parecer — resumo do trabalho autônomo

Branch `correcoes-parecer`, criado a partir do `main` em `df28150`. Cada
tarefa tem uma pasta com script(s), resultados, log e um `RELATORIO.md`
detalhado.

**Regras respeitadas (verificado no fim):** `main` local e remoto continuam
em `df28150`; a tag `v1.0.0` não foi tocada; nenhum arquivo fora de
`correcoes_parecer/` mudou (pacote, `scripts/`, `final_lam02.csv`,
`benchmark_cifsup_lam02.csv` e figuras idênticos ao `main`); nenhum default
foi alterado; packing, χ² e o artigo não foram tocados.

Todo script que usa o método publicado tem uma checagem embutida: sua
versão "atual" reproduz o CSV publicado com 0 diferenças.

## Tarefa 1 — número de inícios (`t1_inicios/`)

**Pergunta:** com quantos inícios o ICP-SAXS deixa de ser
significativamente pior que o cifsup NSD em Chamfer?

**Resposta: com 100.** Semente publicada, Wilcoxon pareado nas 50 entradas:

| Inícios | Chamfer (Å) | p vs cifsup NSD (4,997 Å) |
|---|---|---|
| 3 (publicado) | 5,170 | 9,4e-6 |
| 10 | 5,053 | 3,2e-6 |
| 30 | 5,008 | 0,017 |
| 100 | 4,994 | 0,066 |

* Mesma conclusão nas sementes 1042 e 2042 (100 inícios: p = 0,065 e 0,088).
* **"Não significativo" não é "equivalente":** com 100 inícios o cifsup
  ainda vence em 33 das 50 entradas. Não há correção para comparações
  múltiplas; um teste de equivalência (TOST) seria necessário para afirmar
  mais.
* **Custo:** 100 inícios custam ~33 vezes mais que 3.
* **Achado lateral:** na comparação publicada, em 9 entradas o Chamfer do
  ICP-SAXS usa o envelope amostrado a 5000 pontos e o do cifsup, o envelope
  completo. Corrigir isso muda o nosso Chamfer ~0,005 Å, e as conclusões
  não mudam.

## Tarefa 2 — pesos do Procrustes vs Eq. 1 (`t2_procrustes/`)

**O parecer está certo.** O custo avaliado é a Eq. 1, mas os pesos do
Procrustes não (`src/saxs_icp/core.py`, linhas 129, 140–145 e 161):

* a razão vazamento/preenchimento efetiva é N_P(1+λ)/N_E (mediana 0,65, de
  0,06 a 3,46), não λ = 0,2;
* o vazamento tem um "1 +" fixo;
* o peso cresce com a distância, o oposto do coerente com uma média de
  distâncias.

Comparação num script separado, com a Eq. 1 tirada do README (o .docx não
está nesta máquina):

| Variante | Chamfer, 3 inícios | Chamfer, 30 inícios | custo J, 30 inícios |
|---|---|---|---|
| atual | 5,170 | 5,008 | 6,920 |
| `eq1_ls` (pesos da Eq. 1) | 5,202 (p = 0,28) | 4,987 (p = 0,17) | **6,855 (menor em 42/50, p = 3,8e-7)** |
| `eq1_irls` | 5,406 (pior) | 5,087 (pior) | 6,984 |

Pesos coerentes otimizam melhor o objetivo declarado sem mudar o Chamfer de
forma significativa. Trocar exigiria refazer o benchmark e a calibração de
λ.

## Tarefa 3 — ground truth corrigido (`t3_ground_truth/`)

`ground_truth_v2.py` é uma cópia com as 3 correções:
* (a) o ângulo é calculado só para rotações próprias;
* (b) as sementes são reprodutíveis: crc32 no lugar de `hash()` e `--seed`
  no DAMMIF, verificado byte a byte;
* (c) o cifsup NSD/ICP roda no mesmo teste.

O original está intacto, e os envelopes foram versionados para reuso.

**Resultado desfavorável** (8 entradas × 4 tentativas, sucesso = RMSD < 5 Å):

| Método | Sucesso |
|---|---|
| ICP-SAXS, 3 inícios (publicado) | **1/32** |
| ICP-SAXS, 30 / 100 inícios | 5/32 / 4/32 |
| cifsup NSD / ICP | 8/32 / 12/32 |

* Com ≥ 30 inícios, a pose achada encaixa **melhor** que a "verdadeira" em
  30–32/32: as falhas deixam de ser da busca.
* A "pose verdadeira" deste teste é definida por PCA e é aproximada (RMSD de
  4–10 Å; em 3 das 8 entradas o DAMMIF parece ter produzido a mão trocada).
* **Hipótese não verificada:** se o cifsup inicializa por eixos de inércia,
  como o SUPCOMB, o teste o favorece por construção.

## Tarefa 4 — inicialização por eixos de inércia (`t4_inercia/`)

**Teste sintético:** PCA com 4 inícios recupera 97,7% das poses, contra
36,3% do publicado (McNemar p = 3,8e-54). Funciona em todas as sementes
(94–100%). Mas esse teste é o caso ideal para PCA, porque o envelope é
feito do próprio modelo.

**Dados reais**, que acrescentei por causa desse viés:

| Teste | Publicado | PCA | Referência |
|---|---|---|---|
| Benchmark, Chamfer | 5,165 (3 inícios) | 5,024 com 4 inícios (p = 0,003 vs publicado) ≈ 10–30 inícios aleatórios; 5,000 com 30 | cifsup NSD 4,997; PCA 4 ainda pior (p = 0,001) |
| Envelopes DAMMIF (T3) | 1/32 | 8/32 com 4 inícios | cifsup NSD 8/32 (mesmo viés de PCA da T3) |

PCA é uma forma barata de obter o ganho de muitos inícios, não um ganho
além disso.

## Conclusões para você decidir

1. **A busca é o principal problema do método publicado.** Mais inícios ou
   inicialização por PCA melhoram o Chamfer de forma significativa (5,17 →
   ~5,00 Å) e levam o ICP-SAXS a empatar estatisticamente com o cifsup NSD,
   sem superá-lo.
2. **Os pesos do Procrustes não seguem a Eq. 1.** Corrigir é coerente e não
   piora; o efeito no Chamfer não é significativo.
3. **No ground truth com envelopes do DAMMIF, o cifsup vai melhor.** Parte
   disso pode vir de como a "verdade" é definida (PCA). Esse teste, do
   jeito que está, não serve para uma comparação justa entre os métodos.
4. **Qualquer troca de método** (inícios, PCA, pesos) exige refazer o
   benchmark, a calibração de λ, a comparação com o cifsup e as figuras.

## Pronto para revisar

| Commit | Tarefa |
|---|---|
| `a05e883` | Tarefa 1 |
| `8788f8d` | Tarefa 2 |
| `9f1028a` | Tarefa 3 |
| `a45fc8f` | Tarefa 4 |
| (este) | Resumo |

Para ver o conjunto: `git diff main..correcoes-parecer --stat`.

Nenhuma tarefa travou. Duas coisas que fiz além do pedido, ambas sem mexer
no método:
* a Tarefa 3 ganhou uma opção `--restarts` e rodadas com 30 e 100 inícios;
* a Tarefa 4 ganhou o complemento em dados reais.
