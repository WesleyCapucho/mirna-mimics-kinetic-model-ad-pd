# Modelo cinético de miméticos de microRNA na doença de Alzheimer e na doença de Parkinson: dados, código e resultados

🇬🇧 English version: [README.md](README.md)

Este repositório guarda os dados, scripts, tabelas e figuras de um estudo cinético de miméticos de microRNA (miRNA) para dois eixos: miR-29, BACE1 e Aβ42 na doença de Alzheimer (DA), e miR-7 e α-sinucleína na doença de Parkinson (DP). Não contém texto de manuscrito. Todo número que um manuscrito cita é lido de uma tabela de resultados por um script, e o `scripts/08_verify_consistency.py` recalcula os números derivados a partir das fontes.

## O que contém

- **Uma tabela de parâmetros com origem rastreável** (`data/extracted/kinetic_parameters.csv`): constantes de taxa, meias-vidas e concentrações, cada uma ligada à frase ou à tabela suplementar de onde foi lida, com o sistema em que foi medida.
- **Um modelo EDO por doença**, calibrado só com valores medidos, e uma varredura da única grandeza não medida, a força de repressão, ao longo de seis ordens de grandeza.
- **O custo da dosagem repetida**: o pico que um mimético precisa atingir para manter uma média-alvo, e o que isso faz com o mRNA e a proteína-alvo.
- **Um mapa do espaço de projeto** da dose necessária sobre a força de repressão e o intervalo de dose, para nove séries.
- **Renders de estrutura** de quatro estruturas depositadas (Argonauta 2 com guia e alvo, BACE1, fibrilas de Aβ42 e de α-sinucleína), só para ilustração.

## Relação com o repositório companheiro

O repositório companheiro, [`-Epigen-tica-e-micrornas-em-doen-as-neurodegenerativas`](https://github.com/WesleyCapucho/-Epigen-tica-e-micrornas-em-doen-as-neurodegenerativas) (Zenodo, [10.5281/zenodo.23438256](https://doi.org/10.5281/zenodo.23438256)), guarda uma revisão sistemática e meta-análise de miRNAs circulantes como biomarcadores diagnósticos na DA e na DP. Os dois estudos não compartilham dados nem resultados. Até este repositório existir, os scripts mecanísticos ficavam no companheiro, então as versões arquivadas dele até a 1.2.1 ainda os contêm; daqui em diante este é o lugar deles. A numeração dos scripts continua a do repositório companheiro: o `scripts/12` daqui é o `scripts/12` de lá.

## Estrutura

```
├── data/
│   ├── extracted/kinetic_parameters.csv   # A tabela de parâmetros, uma linha por parâmetro por fonte
│   └── raw/
│       ├── kinetics_2026/                 # Tabelas suplementares de artigos publicados (ver LICENSE-DATA.md)
│       └── structures_2026/               # Estruturas depositadas: 6N4O, 4D8C, 5OQV, 6CU7
├── scripts/
│   ├── 12_ode_models_calibrated.py        # Os dois modelos, calibrados com valores medidos; varreduras dos parâmetros livres
│   ├── 13_structure_figures.py            # Lê as quatro estruturas e registra a proveniência delas
│   ├── 17_mimic_dosing_feasibility.py     # Custo pico/média da dosagem repetida por meia-vida do mimético
│   ├── 36_mimic_target_response.py        # Dosagem, mRNA e proteína juntos; força de repressão; escala de nucleação
│   ├── 40_mimic_design_space.py           # Dose necessária sobre força de repressão e intervalo de dose, nove séries
│   ├── 46_mimic_series_ranking.py         # Qual série mantém um alvo baixo pelo maior intervalo de dose; todo limiar e toda âncora ordenados
│   ├── 47_mimic_series_figure.py          # Figura da ordenação das séries: renders dos alvos acima dos gráficos de intervalo sobre rho (figura do texto principal)
│   ├── 48_mechanism_algorithm_scheme.py   # Esquema do mecanismo molecular das duas doenças e dos seis passos do algoritmo, com renders e gráficos pequenos desenhados das tabelas de resultado
│   ├── 39_structure_renders_article2.py   # Renders do PyMOL das quatro estruturas (só ilustração)
│   ├── 37_manuscript2_numbers.py          # Todo número que um manuscrito cita, lido das tabelas (só números)
│   ├── 38_mechanistic_figures.py          # As sete figuras do texto principal, ajustadas aos limites de arte da revista
│   ├── 18_graphical_abstract.py           # Figura esquemática dos dois eixos, números calculados
│   ├── 19_mechanism_animations.py         # GIFs de mecanismos já simulados
│   ├── 20_structure_story_panel.py        # Compõe as renders de estrutura numa figura só
│   ├── 45_check_figure_artwork.py         # Contraste do texto, espessura de linha, tamanho de página e resolução de imagem em PDFs de figura
│   ├── 08_verify_consistency.py           # Recalcula cada número derivado a partir da fonte; tem de passar antes de versionar
│   └── _bilingual.py, _journal_fit.py, _provenance.py, _viz_mech.py, _viz_style.py   # Módulos de apoio
├── results/
│   ├── tables/                            # Resultados dos modelos, varreduras, espaço de projeto, classes de proveniência, números do manuscrito
│   └── figures/                           # Toda figura em inglês e pt-BR (PNG, PDF, EPS), animações, renders de estrutura
├── docs/                                  # METHODS e DATA_DICTIONARY, em inglês (docs/en) e português (docs/pt-BR)
├── CHANGELOG.md, CHANGELOG.pt-BR.md       # O que mudou entre versões
├── CITATION.cff, LICENSE, LICENSE-DATA.md
└── requirements.txt
```

## Reproduzindo a análise

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python scripts/12_ode_models_calibrated.py
python scripts/17_mimic_dosing_feasibility.py
python scripts/36_mimic_target_response.py        # precisa do 12 e 17; cerca de 3 minutos
python scripts/40_mimic_design_space.py            # precisa do 12, 17 e 36; cerca de 8 minutos
pip install pymol-open-source                     # só para os scripts/13 e 39; instala para o Python a que o pip pertence, então rode esses dois com esse Python
python scripts/46_mimic_series_ranking.py          # precisa do 40
python scripts/13_structure_figures.py
python scripts/39_structure_renders_article2.py   # precisa da tabela de proveniência do 13; vários minutos de ray tracing sem janela
python scripts/37_manuscript2_numbers.py           # precisa do 12, 17, 36, 39 e 40
python scripts/38_mechanistic_figures.py           # precisa do 12, 13, 17, 36, 39 e 40
python scripts/47_mimic_series_figure.py           # precisa do 39, 40 e 46
python scripts/48_mechanism_algorithm_scheme.py     # precisa do 12, 13, 17, 36, 39, 40 e 46
python scripts/19_mechanism_animations.py
python scripts/18_graphical_abstract.py            # precisa do 12, 17 e 13
python scripts/20_structure_story_panel.py
python scripts/08_verify_consistency.py            # tem de passar antes de versionar
```

O `scripts/08` faz cumprir a consistência dos parâmetros, dos resultados dos modelos, da grade do espaço de projeto, das renders e das figuras. Atualmente 991 verificacoes passam.

## Limites

O modelo reafirma o que foi medido e varre o que não foi. Não é uma previsão farmacocinética para nenhum mimético, e as renders ilustram moléculas nos passos dos dois eixos; não testam nada. Ver `docs/pt-BR/METODOS.md`.

## Citação

Ver `CITATION.cff`. O DOI da primeira versão arquivada é acrescentado quando ela for lançada no Zenodo.

## Licenças

Scripts: MIT (`LICENSE`). Dados, tabelas, figuras e documentação: CC BY 4.0 (`LICENSE-DATA.md`), exceto material de terceiros, que mantém os próprios termos; o aviso lista o que é excluído.
