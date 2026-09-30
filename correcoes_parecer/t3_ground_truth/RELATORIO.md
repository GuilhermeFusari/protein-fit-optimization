# Tarefa 3 — Ground truth corrigido

**Resultado principal, desfavorável:** no teste corrigido (8 entradas × 4
tentativas, envelopes do DAMMIF), o ICP-SAXS com os parâmetros publicados
recupera a pose em **1 de 32** tentativas. O cifsup, no mesmo teste,
recupera em 8/32 (NSD) e 12/32 (ICP). Com 30 ou 100 inícios o ICP-SAXS sobe
para 5/32 e 4/32, ainda abaixo do cifsup. Com 30 inícios ou mais, porém, a
pose que ele acha encaixa no envelope **melhor** que a "verdadeira" em 30–32
das 32 tentativas. Pelo critério geométrico que o método otimiza, essas
falhas não são mais falhas da busca.

O `scripts/ground_truth.py` original não foi tocado.

## As três correções (em `ground_truth_v2.py`, cópia do original)

| Defeito | Correção |
|---|---|
| (a) ângulo calculado com reflexões | `angle_proper`: se a transformação total tem det = −1 (o método escolheu a imagem especular), `erro_ang` = NaN e `espelhado` = 1. Essas tentativas saem da média angular, mas continuam contando para RMSD e sucesso. |
| (b) sementes com `hash()` | Semente por entrada = `seed + crc32(código)`. O DAMMIF também recebe `--seed`. Verifiquei rodando duas vezes: envelope do DAMMIF idêntico byte a byte e CSV idêntico (fora a coluna de tempo). |
| (c) cifsup não rodava | O cifsup roda na mesma estrutura embaralhada e no mesmo envelope, métodos NSD e ICP (busca de enantiômeros ligada, o default, como no nosso pipeline). A pose do cifsup é lida do PDB de saída; a matriz aplicada é recuperada por Procrustes com reflexão permitida. |

Nada mais mudou: mesma cadeia CRYSOL → GNOM → DAMMIF, mesma "pose
verdadeira" por PCA, mesmo pipeline próprio (`fit_pipeline`: author,
λ = 0,2, 3 inícios, enantiômeros) e mesmo critério (RMSD < 5 Å). As mesmas
8 entradas e 4 tentativas do teste antigo.

Acrescentei ao v2 duas opções que não mudam o método:

* `--reuse-envelopes`: repete o teste sem rodar o DAMMIF, com os envelopes
  guardados em `envelopes_v2/` (versionados, 764 KB). Os envelopes são
  gravados em PDB com 3 casas decimais, então o reuso difere do cálculo
  original em até 0,07 Å de RMSD. Sucesso e espelhamento saíram idênticos
  nas 32 tentativas.
* `--restarts N`: roda o ICP-SAXS com N inícios nos mesmos envelopes e
  embaralhamentos.

## Resultados

| Método | Sucesso (RMSD < 5 Å) | RMSD mediano (Å) | Espelhado | Erro angular mediano, só rotações próprias | Pose achada encaixa melhor que a verdadeira |
|---|---|---|---|---|---|
| **Antigo** (`ground_truth_final.csv`) | 1/32 | 26,5 (média) | não medido | 93,0° (média, com reflexões) | 26/32 |
| ICP-SAXS, 3 inícios (publicado) | **1/32** | 23,9 | 14 | 31,9° | 18/32 |
| ICP-SAXS, 30 inícios | 5/32 | 14,6 | 9 | 16,5° | 30/32 |
| ICP-SAXS, 100 inícios | 4/32 | 13,8 | 9 | 17,0° | 32/32 |
| cifsup NSD | 8/32 | 10,0 | 12 | 3,8° | 32/32 |
| cifsup ICP | 12/32 | 14,5 | 8 | 5,8° | 28/32 |

Linhas do ICP-SAXS com 3/30/100 inícios: rodadas com envelope reusado
(`ground_truth_v2_reuso_r*.csv`), para comparar os três em condições
idênticas. A linha de 3 inícios dá o mesmo resultado da rodada original
(`ground_truth_v2.csv`: 1/32).

Sucessos por entrada (de 4):

| Entrada | ICP-SAXS 3 | ICP-SAXS 30 | ICP-SAXS 100 | cifsup NSD | cifsup ICP |
|---|---|---|---|---|---|
| SASDUX5 | 1 | 4 | 4 | 4 | 4 |
| SASDWG8 | 0 | 1 | 0 | 4 | 4 |
| SASDUD6 | 0 | 0 | 0 | 0 | 4 |
| outras 5 | 0 | 0 | 0 | 0 | 0 |

## Comparação com o teste antigo

* A taxa de sucesso do ICP-SAXS com 3 inícios é a mesma (1/32), mas os
  envelopes do DAMMIF são outros, porque o original não tinha semente. O
  número bater não significa que as tentativas sejam as mesmas.
* O erro angular antigo (93°, média) misturava reflexões com rotações. Só
  com rotações próprias, a mediana é 31,9°, e 14 das 32 soluções são
  espelhadas.
* O teste antigo não tinha o cifsup, então não dava para saber que o cifsup
  vai melhor neste teste. Agora dá.

## Interpretação e ressalvas

1. **Com 3 inícios, parte das falhas é da busca.** Mais inícios melhoram
   (1 → 5/32; RMSD mediano 23,9 → 14,6 Å), como no teste sintético e na
   Tarefa 1.
2. **Com ≥ 30 inícios, as falhas restantes não são mais da busca.** A pose
   achada encaixa no envelope melhor que a "verdadeira" em 30–32/32. Duas
   explicações possíveis, que este teste não separa:
   * a "pose verdadeira" por PCA é aproximada: o próprio alinhamento por PCA
     deixa RMSD de 4–10 Å entre envelope e estrutura, e em 3 das 8 entradas
     (SASDPT5, SASDV62, SASDVK8) o alinhamento com reflexão encaixa melhor,
     sinal de que o DAMMIF produziu a mão trocada;
   * envelopes de baixa resolução não determinam a orientação: várias poses
     encaixam igualmente bem, e o critério geométrico do ICP-SAXS escolhe
     outra.
3. **Possível viés a favor do cifsup (hipótese, não verificada).** A
   publicação do SUPCOMB, antecessor do cifsup, descreve um alinhamento
   inicial pelos eixos de inércia. Se o cifsup 3.2 faz o mesmo, ele parte
   exatamente do alinhamento que define a "pose verdadeira" neste teste, o
   que o favoreceria por construção. A ajuda do cifsup não diz como ele
   inicializa. Isso vale verificar antes de usar este teste para comparar
   os dois métodos. A Tarefa 4 (inicialização por eixos de inércia no
   ICP-SAXS) é relevante aqui.
4. **Amostra pequena:** 8 entradas × 4 tentativas, com os sucessos
   concentrados em 2–3 entradas. As diferenças entre métodos são
   descritivas; não fiz teste estatístico.

## Arquivos

* `ground_truth_v2.py` — cópia corrigida (o cabeçalho descreve as mudanças).
* `ground_truth_v2.csv` / `.log` — rodada completa: ICP-SAXS (3 inícios) +
  cifsup NSD + cifsup ICP, 96 linhas.
* `ground_truth_v2_reuso_r3.csv`, `_r30.csv` / `.log`, `_r100.csv` / `.log`
  — ICP-SAXS com 3, 30 e 100 inícios nos envelopes guardados.
* `envelopes_v2/` — os 8 envelopes do DAMMIF já alinhados por PCA.

Para reproduzir:

    python3 ground_truth_v2.py --base <pasta com data/sasbdb/> \
        --entries SASDPT5 SASDTS8 SASDUD6 SASDUN8 SASDUX5 SASDV62 SASDVK8 SASDWG8 \
        --trials 4 --out ground_truth_v2.csv --keep-envelopes envelopes_v2
