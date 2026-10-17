# Dicionário de dados

## `data/extracted/kinetic_parameters.csv`

Constantes de velocidade, meias-vidas e concentrações usadas pelo `scripts/12`, uma linha por parâmetro por fonte. 67 linhas de 15 fontes primárias.

| Coluna | Tipo | Descrição |
|---|---|---|
| `param_id` | string | Identificador (K001…). O `scripts/12` carrega os parâmetros por esse id, nunca pela posição |
| `kind` | enum | `numeric` (número impresso na fonte), `qualitative_constraint` (afirmação sem número utilizável, ex.: "indetectável"), `declared_gap` (procurado e não encontrado), `derived` (calculado por nós a partir de uma tabela primária) |
| `axis` | enum | `miR-29/BACE1/Abeta`, `miR-7/SNCA/alpha-synuclein`, `both` |
| `parameter` / `symbol` | string | O que o valor é, e o símbolo usado no modelo |
| `value_as_written` | string | O valor exatamente como impresso na fonte. Em `numeric` e `qualitative_constraint`, precisa aparecer dentro de `verbatim_quote` |
| `value_si` | float | O mesmo valor convertido para as unidades do modelo. Vazio em restrições qualitativas e lacunas |
| `unit_si` | enum | `1/hour`, `M`, `M^-1 s^-1`, `molecules per cell`, `uM`, `um^3`, `fold`, `fractional change`, `fractional increase`, `fraction of patients`, `dimensionless` |
| `unit_as_written` | string | Unidade como impressa na fonte |
| `population` / `condition` | texto livre | Quem ou o que foi medido, e em que condição experimental. Constante de velocidade sem condição é rejeitada pelo `scripts/08` |
| `n` | string | Tamanho amostral como reportado |
| `method` | texto livre | Técnica de medida |
| `species` | string | Organismo ou sistema (LCR humano, linhagem celular, neurônios de rato ou camundongo, proteína recombinante in vitro…) |
| `pmid` / `doi` / `first_author` / `year` / `journal` | | Citação da fonte primária |
| `verbatim_quote` | texto livre | **A frase de onde o valor foi lido.** Vazia de propósito nas linhas `derived` e `declared_gap` |
| `source` | string | Onde o texto foi lido (texto completo no PMC, PDF ou tabela suplementar fornecidos pelo autor deste repositório, texto completo no Scite) |
| `note` | texto livre | Conversões, ressalvas e, nas linhas `derived` e `declared_gap`, como o valor foi calculado ou como a busca foi feita |

**Regras de leitura.**
- Uma linha `declared_gap` não traz valor, citação nem frase. Ela registra que um número foi procurado e não encontrado, e onde. O código das EDOs se recusa a carregá-la.
- Uma linha `derived` **não** é número impresso pela fonte. A derivação está em `note`, e o `scripts/08` a recalcula a partir da tabela arquivada em `data/raw/kinetics_2026/`.
- K038 e K039 registram valores de um modelo que o próprio artigo rejeita. Ficam como aviso e não são usados.
- Um valor impresso numa figura ou na tabela de uma figura é registrado com aquela linha ou rótulo como citação, e o `note` diz de que figura foi lido. Um valor que só existe como posição numa curva não é registrado.
- K052 e K054 são a razão de existir o tipo `derived`. Wilhelm et al. reportam α- e β-sinucleína juntas (K052) e, numa nota de rodapé, a razão entre as duas (K053); nenhuma das duas é, sozinha, uma concentração de α-sinucleína. O K054 combina as duas e diz isso, e o `scripts/08` o recalcula. Citar K052 como α-sinucleína superestimaria o valor em duas vezes.

## `data/raw/kinetics_2026/tushev_2018_table_S1.xls`

Tabela Suplementar S1 de Tushev et al. 2018 (Neuron, DOI 10.1016/j.neuron.2018.03.030): 24.435 isoformas de 3'UTR com símbolo do gene, enriquecimento por tipo celular, localização por compartimento e uma coluna `half.life[hours]`, que são os dados processados para os quais a própria declaração de disponibilidade do artigo aponta. O `scripts/08` procura de novo as linhas de SNCA, BACE1 e APP neste arquivo e recalcula a constante de decaimento agrupada da BACE1. Leia com duas ressalvas, ambas registradas nas linhas que o usam: o decaimento foi observado por 16 h, então um quarto das isoformas traz meias-vidas além da janela e as maiores chegam a 18.858 h; e o resumo do artigo para transcrições de neurônio só se reproduz quando se excluem valores acima de cerca de 25 h, enquanto o da glia não se reproduz de forma alguma.

## `data/raw/kinetics_2026/wilhelm_2014_table_S1.xlsx`

Additional Data Table S1 de Wilhelm et al. 2014 (Science, DOI 10.1126/science.1252884): medidas de immunoblot quantitativo de 64 proteínas presinápticas, com porcentagem da proteína total, número de cópias por sinapse (média ± EPM de quatro preparações), concentração molar e notas de rodapé. A legenda está embutida na planilha como imagem. As concentrações foram calculadas pelos autores sobre o volume sináptico menos o volume mitocondrial; o `scripts/08` confirma isso recalculando o volume implicado por cada par de cópias e concentração e comparando com K050 menos K057, lidos da Figura 1C do mesmo artigo.

## `data/raw/kinetics_2026/schwanhausser_2011_supplementary_table.xls`

Tabela suplementar de Schwanhäusser et al. 2011 (Nature, DOI 10.1038/nature10098): meias-vidas de mRNA e proteína, números de cópias e taxas de síntese em escala genômica em fibroblastos NIH3T3 de camundongo, 5028 linhas. Arquivada sem alteração. Os metadados do arquivo indicam última modificação em novembro de 2012; a errata do artigo (DOI 10.1038/nature11848) saiu em fevereiro de 2013, e não foi confirmado se esta é a versão corrigida. Usada só para as medianas genômicas (K040, K041) e para confirmar a ausência de BACE1 e SNCA (K030).

## `data/raw/structures_2026/`

Coordenadas depositadas, sem alteração: `6N4O.pdb` (Argonauta2 humana com miR-122 e alvo, raio X 2,9 Å), `4D8C.cif` (BACE1 com inibidor sulfona cíclica, raio X 2,07 Å), `6CU7.cif` (fibrila de α-sinucleína completa, polimorfo rod, crio-EM), `5OQV.cif` (fibrila de Aβ(1-42), crio-EM 4,0 Å).

## `results/tables/ode_calibrated_results.json`

Saída do `scripts/12`: os parâmetros livres, suas faixas e o valor central de cada um (medido quando existe medida comparável); as constantes ilustrativas de agregação e por que nenhuma fonte as dá; a comparação entre a dose de mimético exigida e a queda de BACE1 medida; o experimento depuração-versus-produção na DA; o tempo para um mimético de miRNA ser eliminado em cada meia-vida medida; o experimento da comporta de pH da α-sinucleína; as cargas de agregado de Aβ42 humanas expressas em múltiplos de M\*; e um resumo da varredura dos parâmetros livres.

## `results/tables/ode_free_parameter_sensitivity.csv`

Uma linha por ponto da varredura dos parâmetros livres: `free_parameter`, `value`, `abeta_AD_over_control`, `verdict` (`AD_above_control`, `AD_at_or_below_control`, `not_evaluable`). Uma execução que falha numericamente é reportada como `not_evaluable`, nunca contada como inversão. As execuções integram só a cascata do monômero (o bloco ilustrativo de agregação fica de fora, porque não pode mudar o pool de monômero e divergia em dois cantos da varredura); as 18 são avaliáveis.

## `results/tables/structure_figure_provenance.json`

Uma entrada por estrutura, lida do arquivo de coordenadas pelo `scripts/13`: título depositado, método, resolução e critério de resolução, DOI e PMID da citação primária, o papel da figura, e verificações calculadas a partir das coordenadas: para 4D8C, os aspartatos catalíticos encontrados tanto pelo motivo de sequência (DTGS, DSGT) quanto pela distância ao inibidor; para as fibrilas, as cadeias de cada protofilamento, o espaçamento entre camadas e a distância entre protofilamentos. Os números de resíduo são os do depósito.

## `results/tables/mimic_dosing_feasibility.csv`

Uma linha por espécie e intervalo de dose: `species`, `param_id` (a linha da tabela cinética de onde veio a constante de decaimento), `half_life_hours`, `dosing_interval_hours`, `peak_over_average` e `fraction_of_interval_above_half_peak`. Nenhum parâmetro livre e nenhuma dose suposta entram: a dose se cancela na razão entre pico e média.

## `results/tables/mimic_dosing_feasibility.json`

O mesmo cálculo com as comparações que o usam: a constante de decaimento de cada espécie e seu `param_id`, a comparação de dosagem diária contra a referência de estabilidade mediana (penalidade, intervalo equivalente, fator de estabilização necessário) e a suposição declarada de que o mimético entregue é eliminado na taxa endógena — por isso a saída é um fator de estabilização necessário e não um veredito de viabilidade.

## `results/tables/mimic_target_response.csv`

Uma linha por cenário do `scripts/36`: `axis` (AD ou PD), `rho_case` (`scan low`, `scan midpoint`, `scan high`, e para a PD `K065 lower bound`, a repressão mais fraca que o K065 permite para o par), `k_repress`, `rho`, a `species` entregue com seu `param_id` e `half_life_h`, `interval_h`, `protein_case` e `protein_half_life_h`, o `target` e `target_knockdown`, e se é `reachable` com dose até 10^4 vezes o pool basal do repressor. Para um cenário alcançável: a `dose` em bolus e os resultantes `peak_excess`, `peak_pool_over_baseline`, `trough_pool_over_baseline` e `mean_excess` (todos em unidades do pool basal), o `linear_excess` e o `linear_dose` que a forma fechada do estado estacionário pede, `efficiency` (`linear_dose` sobre `dose`, no máximo 1), `knockdown_reached` (a queda média da proteína no ciclo, igual à meta), `protein_ripple` (oscilação pico a vale da proteína sobre sua média), `abeta_ripple` (só AD), `elongation_reduction_m49_8` e `elongation_reduction_m46_0` (só PD: o retardo médio da elongação da fibrila na curva de saturação com K009 e K010, na concentração do botão K054) e `residual` (quão bem o estado periódico fechou). Para um cenário inalcançável as colunas de dose ficam vazias e `knockdown_reached` é o que a maior dose alcança. A dose não depende de `protein_case`: num ciclo periódico a proteína média é igual ao mRNA médio.

## `results/tables/manuscript2_numbers.json`

Gravado pelo `scripts/37_manuscript2_numbers.py`: todo número que o manuscrito do Artigo 2 cita, como mapa plano de chave para texto já arredondado e formatado para o texto, com o arquivo de origem de cada valor em `_sources`. O manuscrito fica fora deste repositório e se refere a estas chaves; uma montagem que não acha uma chave para. Ele seleciona, arredonda e conta a partir das saídas dos `scripts/12`, `17`, `36`, `39` e `40`, do `scripts/_provenance.py` e da tabela de parâmetros cinéticos, e calcula quantidades de aritmética simples (uma razão, uma porcentagem, um tempo até 90%) que rotula como tal. Um assert sustenta cada frase do manuscrito que depende de um padrão nos dados.

## `results/tables/mimic_target_rho_sweep.csv`

O mesmo cálculo num intervalo de dosagem de 24 h sobre a força de repressão ρ de 0,01 a 10⁴ (31 pontos, cinco por década), para todo eixo, espécie entregue e meta: `rho`, se a meta é `reachable`, a `dose`, `peak_pool_over_baseline`, `linear_dose`, `efficiency`, `protein_ripple` e `knockdown_reached`. Existe para que uma afirmação sobre a dose não dependa da varredura declarada de `k_repress`: repressão fraca eleva a dose aproximadamente como 1/ρ, repressão forte se aproxima de um limite, e o menor pico sobre ρ é um limite que vale seja qual for ρ.

## `results/tables/mimic_target_response.json`

O resto do `scripts/36`: as versões de software, os parâmetros de integração, a força de repressão ρ nas pontas e no ponto médio da varredura declarada para cada eixo e, para a DP, no piso do K065 (com o fator de mimético que a forma fechada pede em cada caso), o resumo da varredura de ρ (menor pico diário sobre ρ, e o pico na repressão mais fraca e na mais forte testadas), o escalonamento da taxa de nucleação a partir das ordens de reação medidas do Aβ42, os seis resultados de verificação com seus erros e tolerâncias, as suposições em inglês e português, e as contagens de cenários e de cenários inalcançáveis.

## `results/tables/mimic_design_space.csv` e `mimic_design_space.json`

Gravados pelo `scripts/40_mimic_design_space.py`. O CSV tem uma linha por ponto da grade, 8.379 ao todo (nove séries, 49 valores de ρ, 19 intervalos de dose): `axis` (DA ou DP), `species` com seu `param_id` e `half_life_h`, o `target` e o `target_knockdown`, `rho`, `interval_h`, se a meta é `reachable` com dose até o teto, o `peak_pool_over_baseline` que a alcança e, onde não é alcançável, o `knockdown_reached` no teto. O JSON traz a grade (`rho_grid`, `interval_grid_h`), os três tetos de referência (`thresholds`: 3, 10 e 30 vezes o pool basal), o maior intervalo dentro de cada teto para cada série nos extremos baixo, médio e alto da varredura declarada de ρ e, na DP, no piso do K065 (`interval_within_threshold`, com uma marca: `interp`, `beyond` quando o maior intervalo da grade ainda está dentro do teto, `none` quando nem o menor está), a contagem de `unreached_points` e o bloco `verification` (concordância com o `scripts/36` a 24 h, monotonicidade no intervalo e a conferência do menor intervalo sem alcance contra a forma fechada).

## `results/tables/structure_renders_article2.json`

Gravado pelo `scripts/39_structure_renders_article2.py`: para cada render de PyMOL desenhado para as figuras do estudo mecanístico (em `results/figures/structures/a2_*.png`), a entrada do PDB, o tamanho em pixels do recorte e o que foi selecionado (cadeias, nucleotídeos da semente, aspartatos catalíticos). Para as vistas laterais das duas fibrilas, `layer_transformation` lista, por protofilamento, a cadeia que foi movida, a cadeia a que foi ajustada, o número de pares de C-alfa, o rms do ajuste, a subida por camada (angstrom), o giro por camada (graus) e as cadeias ao longo do eixo. A camada translúcida depois da última cadeia é construída a partir dessa transformação; é ilustração, não resultado.

## `results/tables/parameter_provenance_classes.csv`

Gravado pelo `scripts/_provenance.py` (por meio dos `scripts/37` e `38`): uma linha por linha da tabela de parâmetros cinéticos, com `param_id`, `kind`, `axis` (AD, PD ou both), `source` (primeiro autor e ano, vazio para a lacuna declarada), `system` (onde foi medida: `people`, `human_pm`, `animal_brain`, `neurons`, `cells`, `cell_free` ou `none`) e `block` (o que alimenta: `mirna`, `mrna`, `repression`, `protein`, `flux` ou `aggregation`). As classes são declaradas no módulo e conferidas contra espécie, população, condição e método de cada linha; o `scripts/08` exige que a tabela gravada seja igual à classificação.

## `results/tables/graphical_abstract_values.json`

Todo número desenhado em `results/figures/graphical_abstract.*.png`, carregado do `ode_calibrated_results.json`, do `mimic_dosing_feasibility.json` e do `kinetic_parameters.csv` no momento do desenho, e não digitado como string separada, mais os `param_id`s em que cada afirmação da figura se apoia. Nada aqui é cálculo novo; é o registro de qual número já verificado foi parar onde.

## `results/tables/mechanism_animations.json`

Os números que cada animação em `results/figures/anim_*.gif` desenha: a razão final Aβ42 DA/controle, os tempos de eliminação das duas espécies, e as razões pico sobre média das duas espécies e a penalidade entre elas. Todo campo é recalculado pelo `scripts/08` a partir do `ode_calibrated_results.json` e do `mimic_dosing_feasibility.json` e exigido a bater, para que o manifesto não possa se afastar das tabelas de onde seus próprios números vieram.

## `results/figures/anim_*.gif`

Três GIFs animados, cada um nos dois idiomas: `anim_ad_monomer` (Aβ42 acumulando até o estado estacionário, sem parâmetro livre), `anim_mimic_washout` (um bolus único de mimético decaindo na meia-vida medida), `anim_dosing_sawtooth` (a penalidade de pico sobre média da dosagem repetida como um ciclo em movimento). Cada um é produzido pelo `scripts/19_mechanism_animations.py`, que importa suas equações do `scripts/12` e do `scripts/17` em vez de rederivá-las. A trajetória da comporta de pH da α-sinucleína não está entre eles, porque sua constante de velocidade em pH ácido é ilustrativa e não medida (K042).

## `results/figures/`

Toda figura existe duas vezes, `<nome>.en.png` e `<nome>.pt-BR.png`. As duas são desenhadas pelo mesmo código, a partir dos mesmos vetores, na mesma execução; só o texto dos rótulos muda.

## `results/figures/structure_story_panel.*.png`

Os quatro renders do `scripts/13` compostos numa figura só, com os fatos que acompanham (PDB id, método, resolução, resíduos da díade catalítica, número e espaçamento dos protofilamentos) tirados do `structure_figure_provenance.json`, não redigitados. Nenhum render novo, nenhuma afirmação nova - só o layout e o agrupamento são novos.

## `results/tables/mimic_series_ranking.csv` e `mimic_series_ranking.json`

Gravados pelo `scripts/46`. Uma linha por série, âncora e limiar: `axis`, `target`, `target_knockdown`, `species`, `param_id`, `half_life_h`, `anchor` (scan low, midpoint, high ou K065 floor), `rho`, `threshold_x` (3, 10 ou 30), `interval_h` (maior intervalo de dose dentro do limiar; 0 quando nenhum serve, 168 quando todos os da grade servem), `flag` (`interp`, `none` ou `beyond`), `reference` (True para a referência de miRNA mediano) e `rank_among_studied` (ordenação por competição entre espécies medidas; vazio para a referência). O JSON declara o critério nos dois idiomas, a série em destaque de cada eixo, se o primeiro lugar muda e os empates.

## `data/raw/sequences_2026/` e `data/extracted/utr3_sequences.fa`

Entradas do `scripts/49`, com a proveniência em `provenance.json` (URL, `Last-Modified` quando o servidor dá, tamanho, SHA-256 de cada arquivo e sequência). `mirbase_mature_selected.fa` tem as cinco sequências maduras usadas (por acesso: MIMAT0000086, MIMAT0000100, MIMAT0000681, MIMAT0000252, MIMAT0000439) como escritas no `mature.fa` do miRBase (o arquivo não diz a versão; o nome do MIMAT0000439 ali é hsa-miR-153-1-3p, chamado miR-153-3p nas tabelas). `NM_012104.6.gb` e `NM_000345.4.gb` são os dois registros RefSeq (ambos MANE Select) como arquivos GenBank do NCBI. `utr3_sequences.fa` tem o 3′UTR de cada registro (a sequência depois do fim da feature CDS, de 5′ a 3′ do mRNA, escrita com U; o cabeçalho dá o trecho no registro). `refGene_selected.tsv` tem as duas linhas da tabela refGene do UCSC para hg38, usadas só numa conferência: o 3′UTR do SNCA lido do hg38 é idêntico, e o do BACE1 é o 3′UTR do RefSeq mais 12 nt na ponta 3′ (resultado em `provenance.json`).

## `results/tables/mirna_site_pairing.csv` e `.json`

Escritos pelo `scripts/49`. Uma linha por sítio de pareamento da semente de um microRNA num 3′UTR: `gene`, `transcript`, `utr3_length_nt`, `mirna`, `mirna_accession` (MIMAT), `mirna_length_nt`, `site_type` (8mer, 7mer-m8 ou 7mer-A1), `utr_start` e `utr_end` (contados de 1 no primeiro nucleotídeo depois do códon de parada), `a1_position` (o nucleotídeo do 3′UTR que fica diante do nucleotídeo 1 do microRNA), `seed_pairs_2_8`, `supplementary_pairs_13_16`, `pairs_total` e `wobble_pairs` (sem lacunas, Watson–Crick mais G–U), `target_5to3` e `pairing_marks` (`|` Watson–Crick, `:` G–U). É uma contagem de complementaridade, não uma energia de ligação e não prova de que o sítio é funcional.
