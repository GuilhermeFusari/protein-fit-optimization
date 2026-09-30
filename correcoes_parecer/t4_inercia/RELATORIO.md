# Tarefa 4 — Inicialização por eixos de inércia

**Resposta curta:** melhora, e muito, no teste sintético: com 4 inícios por
eixos de inércia a pose é recuperada em **97,7%** dos casos, contra 36,3%
com a inicialização publicada (3 inícios) e 38,7% com o mesmo orçamento de
4 (McNemar exato p = 4,7e-52). Esse teste, porém, é o caso ideal para a
inicialização por PCA, porque o envelope é construído a partir do próprio
modelo. Nos **dados reais** o ganho é menor: no benchmark das 50 entradas,
PCA com 4 inícios equivale a 10–30 inícios aleatórios em Chamfer. Não passa
do platô de ~5,00 Å que muitos inícios aleatórios já atingem, e continua
significativamente pior que o cifsup NSD.

Nada no pacote foi alterado.

## O que foi implementado (`inicializacao_inercia.py`)

Cópia do laço de `core.align` com uma única mudança: a escolha das
orientações iniciais de cada mão (original e espelhada).

* **atual** (publicado): a pose de entrada + (k − 1) rotações aleatórias.
* **pca**: os eixos principais do modelo alinhados aos do envelope, nas 4
  combinações de sinal com det = +1 e centroide sobre o do envelope, mais
  (k − 4) rotações aleatórias.

O ICP, o custo, a escolha da melhor pose e a busca de enantiômeros são os
do pacote. Protocolo idêntico ao de `scripts/teste_busca.py` (mesmos
envelopes, embaralhamentos e sementes). **Checagem:** a variante `atual`
com 3/10/30 inícios reproduz as 900 linhas de `benchmark/teste_busca.csv`
com 0 diferenças.

## Teste sintético (pose verdadeira exata; 50 entradas × 6 sementes)

| Inicialização | Inícios por mão | Acerto | Falhas com custo > verdadeiro |
|---|---|---|---|
| atual (publicado) | 3 | 36,3% (109/300) | 182/191 |
| atual | 4 | 38,7% (116/300) | 171/184 |
| atual | 10 | 72,7% (218/300) | 75/82 |
| atual | 30 | 94,7% (284/300) | 11/16 |
| **pca** | **4** | **97,7% (293/300)** | 7/7 |
| pca | 10 | 97,7% (293/300) | 6/7 |
| pca | 30 | 97,3% (292/300) | 5/8 |

Comparações pareadas (McNemar exato, mesmas 300 tentativas):

| Comparação | Só pca acerta | Só atual acerta | p |
|---|---|---|---|
| atual 4 vs pca 4 (mesmo orçamento) | 178 | 1 | 4,7e-52 |
| atual 3 (publicado) vs pca 4 | 185 | 1 | 3,8e-54 |
| atual 10 vs pca 10 | 76 | 1 | 1,0e-21 |
| atual 30 vs pca 30 | 8 | 0 | 0,008 |

Por semente, pca com 4 inícios fica entre 94% e 100%. A inicialização
atual com 3 inícios varia de 8% a 80%: a dependência do sorteio praticamente
desaparece. Nenhuma falha com pca escolhe a mão errada. As 7–8 falhas
restantes se concentram em 7 modelos (SASDTS8, SASDUZ5, SASDV75, SASDV85,
SASDV95, SASDVG6, SASDX24), em parte os mesmos casos de forma ambígua do
`teste_busca.py`.

**Por que este teste favorece a PCA:** o envelope sintético é o próprio
modelo engrossado, então os eixos principais dos dois coincidem quase
exatamente, e um dos 4 inícios por PCA já cai perto da pose verdadeira. Em
envelopes reais (reconstruções ab initio de baixa resolução) os eixos só
coincidem aproximadamente. Por isso rodei o complemento abaixo.

## Complemento: dados reais (`inercia_dados_reais.py`)

**A. Benchmark das 50 entradas** (Chamfer contra o envelope completo, como
o cifsup; inicialização atual lida da Tarefa 1, mesma semente 42)

| Inicialização | Inícios | Chamfer (Å) | p vs cifsup NSD | p vs publicado (atual 3) | Melhor que atual 3 |
|---|---|---|---|---|---|
| atual (publicado) | 3 | 5,165 | 1,2e-5 | — | — |
| atual | 10 | 5,048 | 6,4e-6 | | |
| atual | 30 | 5,003 | 0,033 | | |
| atual | 100 | 4,989 | 0,11 | | |
| pca | 4 | 5,024 | 0,0012 | 0,003 | 31/50 |
| pca | 10 | 5,009 | 0,0015 | 5,2e-5 | 34/50 |
| pca | 30 | 5,000 | 0,046 | 9,5e-8 | 40/50 |
| cifsup NSD | | 4,997 | | | |

A PCA com 4 inícios dá quase o mesmo ganho que 30 inícios aleatórios, com
~7 vezes menos cálculo. Mas não ultrapassa o platô: com 30 inícios, PCA
(5,000) e inicialização atual (5,003) empatam na prática, e as duas
continuam no limite da significância contra o cifsup NSD.

**B. Envelopes do DAMMIF da Tarefa 3** (8 entradas × 4 tentativas, mesmos
embaralhamentos)

| Inicialização | Inícios | Sucesso | RMSD mediano (Å) |
|---|---|---|---|
| atual (publicado) | 3 | 1/32 | 23,9 |
| atual | 30 | 5/32 | 14,6 |
| pca | 4 | 8/32 | 9,9 |
| pca | 30 | 6/32 | 13,3 |
| cifsup NSD (Tarefa 3) | | 8/32 | 10,0 |
| cifsup ICP (Tarefa 3) | | 12/32 | 14,5 |

Com PCA, o ICP-SAXS alcança o cifsup NSD neste teste (8/32). **Mas a "pose
verdadeira" deste teste é definida por alinhamento de PCA**, então uma
inicialização por PCA é favorecida por construção. O mesmo vale para o
cifsup, se ele inicializa por eixos de inércia (hipótese da Tarefa 3). Este
resultado é coerente com a hipótese de viés, mas não a prova.

## Interpretação

1. A inicialização por eixos de inércia resolve quase todo o problema de
   busca quando os eixos do envelope refletem os do modelo, e torna o
   resultado praticamente independente da semente.
2. Em envelopes reais ela é uma forma barata de obter o ganho de muitos
   inícios aleatórios. Não é um ganho além disso.
3. Nenhuma das duas mudanças (PCA ou mais inícios) coloca o ICP-SAXS
   claramente à frente do cifsup NSD em Chamfer no benchmark. No máximo,
   o deixa estatisticamente indistinguível (Tarefa 1).
4. A inicialização por PCA com 4 inícios **mais** inícios aleatórios (o
   "eixos de inércia + 30–100 inícios" planejado no KNOWN_ISSUES) é
   defensável: no sintético não perde para nenhuma alternativa, e no real
   iguala o melhor que os inícios aleatórios conseguem. A decisão de
   trocar o método é sua.

## Arquivos

* `inicializacao_inercia.py` → `inicializacao_inercia.csv` (2100 linhas),
  `.log` com a checagem. Leva ~60 min em 8 núcleos.
* `inercia_dados_reais.py` → `inercia_dados_reais.csv`, `.log` (partes A e
  B). Leva ~10 min.
