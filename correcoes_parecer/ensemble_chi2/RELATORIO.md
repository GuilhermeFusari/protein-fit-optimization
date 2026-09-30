# Frente 2 — Validação do packing por χ² contra I(q) experimental

**Resumo:** o pipeline funciona e foi testado no SASDFQ9 (RNA). O ensemble
de 20 membros do packing salvo reproduz a curva experimental com
**χ² = 1,84** (1,11 com constante aditiva). Mas esse packing escolheu **20
de 20** candidatos, que são o ensemble que o EOM depositou já ajustado à
curva. O número valida o cálculo, não a seleção. No teste que valida a
seleção (escolher k < 20), **o packing atual seleciona conformações que
reproduzem a curva pior que subconjuntos aleatórios**: percentil mediano
90–95% da distribuição do acaso, em 10 sementes. A causa aparente é uma
preferência por conformações estendidas.

**GRB2: não foi possível.** A curva experimental da amostra do envelope
(GRB2 W60, pH 7) não está em nenhum lugar do disco. Há curvas de outras
amostras da GRB2 (W121A, Y160F), mas usá-las seria comparar com outra
proteína; não improvisei.

Nada do packing publicado foi alterado: `resultado_packing_SASDFQ9/` foi só
lido, e os novos packings rodaram em pastas próprias.

## Dados e ferramentas

| Item | SASDFQ9 | GRB2 (W60 pH 7) |
|---|---|---|
| I(q) experimental | **sim**: `dados/SASDFQ9_experimental.dat` (cópia de `SAXS_shape_prediction/saxs_v2/data/raw/complete/SASDFQ9.dat`), 261 pontos, q = 0,010–0,270 Å⁻¹ | **não encontrado** |
| Envelope | `dados/SASDFQ9_envelope.cif` (cópia de `data/3 TESTES/SASDFQ9/SASDFQ9/pddf/SASDFQ9-1.cif`, o indicado em `LEIAME_DADOS.md`) | `data/envelope.cif` |
| Candidatos | 20 modelos `fit1_model1..20` depositados | 20 mil quadros (`data/pdbs_w60`) |
| CRYSOL | 3.2.1, disponível | — |

O `SASDFQ9_fit1.fit` depositado é um ajuste **EOM/GAJOE** ("ENSEMBLES: 50,
CURVES: 20", "Chi: 1,042"): os 20 modelos já foram selecionados contra essa
curva.

## Pipeline (`ensemble_chi2.py`)

1. **Perfil de cada membro:** CRYSOL 3.2, parâmetros de hidratação default
   (não ajustados), mesmos para todos os membros.
2. **Perfil do ensemble:** média simples dos perfis (pesos iguais). Em
   solução diluída as moléculas espalham de forma independente, então as
   posições das cópias no envelope **não entram** no χ². Rodar o CRYSOL no
   arquivo combinado `ALL_TOP_ALIGNED.pdb` estaria errado, porque trataria
   as cópias como um único complexo.
3. **χ² reduzido** contra o I(q) experimental, com fator de escala ajustado
   analiticamente, e uma variante com constante aditiva.
4. **Controle:** o χ² de cada membro isolado.

**Ajuste de formato necessário:** os modelos do SASDFQ9 usam o nome de
resíduo `rU` e nomes de átomo antigos/ambíguos (`O1P`, açúcar sem apóstrofo:
`C1` em vez de `C1'`), e o CRYSOL 3.2 abortava. O script corrige numa cópia
temporária (`rU` → `U`, ocupância vazia → 1,00) e roda com
`--explicit-hydrogens=yes`, em que o CRYSOL usa só o elemento de cada átomo.
Os 463 átomos de cada modelo são lidos (0 desconhecidos, conferido no log).
Consequência: os hidrogênios não entram (os arquivos não os têm), o que é
um erro sistemático igual para todos os membros.

## Resultado 1: o ensemble do packing salvo (20 de 20)

![perfil](resultado_SASDFQ9/perfil_ensemble_vs_experimental.png)

| | χ² (só escala) | χ² (escala + constante) |
|---|---|---|
| **Ensemble dos 20** | **1,84** (χ = 1,35) | 1,11 |
| Melhor membro isolado (model16) | 4,69 | |
| Mediana dos membros isolados | 14,6 | |
| Pior membro isolado (model18) | 64,8 | |
| Referência: ajuste EOM depositado | "Chi" = 1,04 | |

O ensemble reproduz a curva muito melhor que qualquer membro isolado, e
chega perto do valor do EOM depositado. A diferença é esperada, porque o EOM
ajusta a hidratação e usa pesos. Isso **valida o pipeline**: com o ensemble
certo, o χ² sai próximo de 1. **Não valida o packing**, porque o packing
salvo pegou todos os 20 candidatos.

## Resultado 2: a seleção do packing contra o acaso

Teste (`selecao_vs_aleatorio.py`): do pool de 20, o `saxs-icp pack` atual
(λ = 0,2) escolhe k = 5 ou k = 10 cópias, com 10 sementes. O χ² da seleção
é comparado com o de **todos** os subconjuntos de k membros (15.504 para
k = 5, 184.756 para k = 10). A seleção usa só a geometria; a curva
experimental não entra nela, então o χ² é uma validação independente.

| k | χ² da seleção do packing (mediana de 10 sementes) | χ² do acaso: mediana (p5–p95) | Melhor subconjunto possível | Percentil do packing (mediana; mín–máx) |
|---|---|---|---|---|
| 5 | **9,05** | 3,80 (1,62–10,92) | 1,13 | **90,5%** (71,9–98,5) |
| 10 | **6,21** | 2,20 (1,31–6,26) | 1,12 | **94,8%** (58,5–98,3) |

Percentil = % dos subconjuntos aleatórios com χ² menor. Nas 20 execuções,
a seleção do packing foi pior que a mediana do acaso em todas.

**Por quê (diagnóstico):** o packing prefere conformações estendidas. O Rg
médio dos escolhidos (k = 5) é 31,0 Å, contra 27,6 Å do pool e 26,9 Å do
envelope. Os modelos 15, 12 e 13 (Rg 36–38 Å, χ² isolado 29–33) são
escolhidos com frequência, e os compactos de melhor χ² quase nunca. Uma
explicação plausível, não testada: o termo de preenchimento do custo
recompensa estruturas que alcançam as extremidades do envelope, e o
envelope é uma forma média, enquanto a curva é dominada pelas formas mais
compactas.

## Ressalvas

1. **Um único sistema, pool pequeno e já pré-selecionado.** Os 20 candidatos
   já são uma seleção EOM feita para ajustar a curva. Não sei se o resultado
   se repete com um pool grande e não filtrado, como os quadros de dinâmica
   da GRB2, para os quais falta a curva.
2. **Hidratação fixa e sem hidrogênios.** Os χ² absolutos dependem disso. A
   comparação entre subconjuntos usa os mesmos perfis, então é menos
   afetada.
3. **Pesos iguais.** O packing não dá pesos às cópias. O EOM usa
   multiplicidades.
4. **O packing testado é o código atual** (λ = 0,2), não o que gerou os
   números publicados (versão anterior, ver `KNOWN_ISSUES.md`).

**Leitura:** isso não decide nada sobre o packing publicado. A decisão
científica é sua e do seu orientador. Mas no único caso que pude testar, o
critério geométrico do packing não seleciona conformações compatíveis com o
I(q) melhor que o acaso.

## O que falta para testar a GRB2

A curva I(q) experimental da amostra GRB2 W60 pH 7, a mesma usada para
gerar o envelope `data/envelope.cif` (GASBOR, "W60_pH7"). Com ela, o mesmo
pipeline roda no packing de 30 cópias.

## Arquivos

* `ensemble_chi2.py` — pipeline (ensemble → CRYSOL → χ² + gráfico).
* `selecao_vs_aleatorio.py` — seleção do packing contra todos os
  subconjuntos.
* `dados/` — cópias da curva experimental e do envelope do SASDFQ9.
* `pool_SASDFQ9/` — os 20 candidatos (cópias dos `RANK_*.pdb` do packing
  salvo, sem o prefixo; a posição não afeta o espalhamento).
* `resultado_SASDFQ9/` — `resultado.json`, `chi2_membros.csv`,
  `perfis_membros.dat`, `perfil_ensemble.dat`, gráfico, log.
* `selecao_SASDFQ9/` — `selecao_vs_aleatorio.json` e o `report.txt` de cada
  uma das 20 execuções do pack. Os PDBs alinhados não foram versionados,
  porque são regenerados com as sementes.
