# Frente 1 — Existem entradas com pose verdadeira real?

Só investigação: nenhum teste de recuperação foi rodado. Aguarda o seu OK.

**Resposta curta:** existem **95 entradas** em que o modelo atômico e o
envelope depositados já estão no mesmo sistema de coordenadas: 7 no
benchmark dos 50 e 88 no `Dataset_Triado/3_Fila_Dammin`. Mas **nenhuma tem
uma superposição comprovadamente independente de cifsup/SUPCOMB ou PCA**:

* 3 das 7 do benchmark têm no próprio arquivo o registro de que foram
  superpostas com **SUPCOMB** (SASDUD6), **CIFSUP** (SASDX24) ou
  **ALPRAXIN** (SASDUZ5, alinhamento por eixos principais, ou seja, PCA).
  Pela sua regra, essas saem.
* As outras 92 não têm registro de como foram superpostas. A análise abaixo
  mostra que **alguém moveu o modelo para o referencial do envelope** (não é
  coincidência nem PCA independente), mas o método não está nos arquivos
  nem na documentação do SASBDB que consultei. A ferramenta padrão para isso
  no ATSAS é o SUPCOMB/cifsup, então é plausível que parte ou todas essas
  poses venham dele. **Não consigo confirmar nem excluir.**

A decisão de aceitar essas poses como verdade é sua.

## Como verifiquei (`investigar_referencial.py`, `triagem_dataset_triado.py`)

Para cada par (modelo, envelope), **na pose depositada, sem mover nada**:

1. **Centroides:** distância entre o centroide dos Cα e o do envelope.
2. **Orientação:** o Chamfer da pose depositada comparado com o de rotações
   aleatórias do modelo em torno do próprio centroide (1000 por par no
   benchmark, 200 na triagem). Centroides coincidentes sozinhos não provam
   nada, porque envelopes do DAMMIF e muitos modelos vêm centrados na
   origem. Se a orientação depositada for melhor que ≥ 99% das aleatórias,
   ela é especial: já encaixa.
3. **Registro de superposição:** busca por SUPCOMB, CIFSUP, ALPRAXIN,
   DAMSUP, SUPALM e "superimposed" no texto do envelope.
4. **PCA disfarçado** (`diagnostico_eixos_principais.csv`): se modelo e
   envelope estivessem **cada um** nos próprios eixos principais, a
   coincidência viria de dois alinhamentos por PCA independentes. Testei
   pela correlação entre os eixos x, y, z das coordenadas de cada objeto
   (≈ 0 quando o objeto está nos próprios eixos principais; em 2000
   orientações aleatórias isso nunca ocorreu).

Critério "mesmo referencial": centroides a < 5 Å **e** orientação
depositada melhor que ≥ 99% das aleatórias.

## Resultados

| Conjunto | Pares examinados | Entradas no mesmo referencial | Com registro de SUPCOMB/CIFSUP/ALPRAXIN | Sem registro |
|---|---|---|---|---|
| Benchmark dos 50 | 50 | **7** (+1 limítrofe, SASDXF5) | 3 (SASDUD6, SASDX24, SASDUZ5) | 4 (SASDXF6, SASDWT6, SASDVG6, SASDWS6) |
| `Dataset_Triado/3_Fila_Dammin` | 5271 pares (2347 entradas legíveis) | **88** | 0 | 88 |

A lista completa das 88, com as métricas, está em
`triagem_dataset_triado.csv` (filtro: `sep_centroide < 5` e
`pct_orientacao <= 1`, melhor par por entrada). As 7 do benchmark estão em
`referencial_benchmark50.csv`.

**Nenhuma das 95 é PCA disfarçado:** em nenhuma o modelo e o envelope estão,
ambos, nos próprios eixos principais. O padrão das 88 do `Dataset_Triado` é
outro: o envelope do DAMMIF está na origem com orientação própria (87/88), e
o modelo foi colocado nesse referencial (48/88 com centroide a < 0,5 Å da
origem). Isso é uma superposição deliberada do modelo sobre o envelope.

**Mão:** em várias das 88, a imagem especular da pose depositada também
encaixa melhor que 99% das rotações aleatórias (`pct_espelho` baixo). O
envelope não distingue bem a mão nesses casos. O teste precisaria contar
como acerto a pose depositada ou a sua imagem especular, ou separar os dois
casos.

## Como seria o teste (se você aprovar a verdade)

1. **Verdade:** a pose depositada do modelo (Cα) no referencial do envelope
   depositado do SASBDB, os dois lidos como estão, sem nenhum alinhamento
   nosso.
2. **Embaralhamento conhecido:** rotação uniforme + translação N(0, 15 Å)
   aplicadas ao modelo, com sementes fixas (6 por entrada, como em
   `teste_busca.py`).
3. **Ajuste:** `saxs_icp.align` com os parâmetros publicados (e, à parte,
   com 30 inícios e com a inicialização por PCA da Tarefa 4, só para
   medir).
4. **Métrica:** RMSD pareado dos Cα contra a pose depositada; sucesso se
   < 5 Å. Separadamente: acerto na mão depositada, acerto na mão espelhada,
   e se a pose achada tem custo menor que a depositada (falha da busca ou
   da função).
5. Nada de cifsup, PCA ou DAMMIF para definir a verdade.

## O que você precisa decidir

* **Opção 1:** aceitar as 92 superposições sem registro como verdade
  "depositada", deixando explícito no texto que o método de superposição é
  desconhecido e provavelmente SUPCOMB/cifsup em parte dos casos. É real no
  sentido de não ter sido feita por nós, mas não é independente de
  ferramentas de superposição.
* **Opção 2:** considerar que não há verdade independente no SASBDB e
  definir juntos um plano B.

**Observação lateral:** nas 7 entradas do benchmark no mesmo referencial, o
primeiro início do alinhador (a pose de entrada) já é a pose superposta. Nas
entradas SASDUD6, SASDXF6, SASDWS6 e SASDVG6, o Chamfer publicado é
praticamente igual ao da pose depositada. Isso não invalida o benchmark,
mas nessas entradas o método começou da resposta.

## Arquivos

* `investigar_referencial.py` → `referencial_benchmark50.csv`
* `triagem_dataset_triado.py` → `triagem_dataset_triado.csv`, `.log`
* `diagnostico_eixos_principais.csv` — teste do PCA disfarçado nas 95
  candidatas.
