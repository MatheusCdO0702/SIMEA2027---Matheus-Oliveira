# V22.3 — roteiro de execução e avaliação

Preparada em 17/11/2026 a partir da V22.2 enviada ao professor. **Não foi executada pelo editor.** Foi feita revisão estática; a validação numérica começa com sua rodagem. As versões e os resultados anteriores permanecem guardados.

## 1. Preparar a pasta

Extraia o ZIP inteiro em uma pasta nova, por exemplo `C:\Users\wcorr\Downloads\Simulador_Auto2_V22_3`. Não misture com a V22.2. O pacote contém o cache NASA utilizado anteriormente, com manifesto e hash; não precisa baixar novamente nem renomear seus dois arquivos. Não contém resultados novos.

Abra o PowerShell nessa pasta. Ambiente de referência: Python 3.13.7. Se já usa esse Python com as dependências abaixo, pode manter o ambiente. Para isolar:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Os comandos seguintes usam esse ambiente, sem necessidade de ativá-lo. No ambiente existente, substitua `.\.venv\Scripts\python.exe` por `python`. Esta versão não utiliza pvlib nem executável CBC externo.

## 2. Verificações

```powershell
.\.venv\Scripts\python.exe .\Simulador_V22_3.py --modo verificacao
```

Resolve A–F, G (arbitragem controlada), tolerância de ultrapassagem e conferências financeiras. Referências históricas não são resultados da nova matriz. A extensão de G para 764 kWh/dia é ordem de grandeza, com perdas e 22 dias de ponta; não deve igualar toda a economia C2−C3.

Abra a pasta indicada ao final. Em `tabelas/verificacao.csv`, todas as linhas devem ter `aprovado=True`. Se houver falha ou exceção, **pare nessa etapa** e envie a pasta. Não altere os valores esperados para fazer o teste passar.

## 3. Conferência dirigida

```powershell
.\.venv\Scripts\python.exe .\Simulador_V22_3.py --modo conferencia
```

Inclui verificações; dez cenários-base Verde/Azul; crédito zero com 6/8 h; dois carregadores; dia observado próximo do P10; e financeiro anual C1/C2/C3. Não é a matriz final E5. Só gera figuras dos blocos presentes, não as nove da matriz completa.

Conferir antes de prosseguir:

- Auditoria sem excesso acima do alvo, erro de balanço, violação de janela/taper ou ciclo aberto. Tolerâncias no registro.
- `VIAVEL_COM_GAP` preserva solução disponível sem certificado dentro do gap desejado; não é infactibilidade. `SEM_SOLUCAO` também não comprova infactibilidade. `INVIAVEL_COMPROVADO` exige status do problema estrito; folga diagnóstica não entra na economia.
- Filtro `bloco=base` em `cenarios.csv`: C0, C0b, C1, C2, C3 nas duas modalidades. BESS ocioso com 10 h pode ser legítimo.
- `frota_prefixo_54_100=True` no registro; mesma frota nas comparações. P10 observado tem data registrada.
- Custos reconciliados com objetivo e regularizador. Diferenças pequenas C2−C3 devem ser lidas junto com o gap.

Envie essa rodada para avaliarmos antes da matriz maior, sobretudo se N=2 permanecer com gap elevado. Não aumente todos os tempos sem necessidade.

## 4. Matriz E5

```powershell
.\.venv\Scripts\python.exe .\Simulador_V22_3.py --modo e5
```

Inclui base e frota nas duas modalidades; carregadores, permanência, crédito e perfis; P5–P15 e P10; contrato 400–800/50 kW; financeiro anual; fronteira 4×3 com confirmações MILP; BESS 500/750/1.000 kWh a 4 h; robustez do regularizador; três sementes LP/MILP. Sensibilidades complementares usam Verde; não são todas duplicadas em Azul.

Tempo total desconhecido até testar. Solver começa em 60 s, escalona tentativas e limita o caso estrito a aproximadamente 900 s acumulados, incluindo preparação. Diagnóstico pode acrescentar até 600 s, além da construção do modelo. Mantenha o computador ligado. Não abra CSVs no Excel durante a gravação.

Menu de parâmetros agrupados e edição de controles de execução:

```powershell
.\.venv\Scripts\python.exe .\Simulador_V22_3.py
```

Cada chamada cria pasta com data/hora/modo. Não retoma rodadas antigas, evitando mistura de códigos/parâmetros. Se interromper, blocos concluídos ficam identificados como **parciais**, não como E5 completa.

## 5. Organização

```text
resultados/E5_DATA_HORA_MODO/
  LEIA_PRIMEIRO.md
  entrada/
    frota_base.csv
    perfis_solares.csv
  tabelas/
    verificacao.csv
    cenarios.csv
    comparacoes.csv
    financeiro_anual.csv
    fronteira.csv
    capacidade_bess.csv
    perfis_operacao.csv
  figuras/
    01_operacao_base.png
    02_frota.png
    03_carregadores.png
    04_permanencia.png
    05_fronteira_anual.png
    06_contrato.png
    07_sementes.png
    08_baixa_geracao.png
    09_capacidade_BESS.png
  auditoria/
    auditoria.csv
    solver.csv
    registro_execucao.json
    codigo_executado.py
```

Só são criados arquivos de blocos executados. Cache solar integral fica uma vez, fora das rodadas. Não há log bruto por solver, NPZ por caso nem figuras duplicadas em PNG/PDF. Os dois CSVs de auditoria comprovam conservação e status/gap/tempo e precisam ser guardados.

CSV: separador `;`, decimal `,`, UTF-8; `N/D` significa indisponível, nunca zero. Filtre `cenarios.csv` por `bloco`. LP relaxado não certifica agenda inteira de carregadores.

## 6. Avaliação completa e retorno

1. Conferir conclusão/verificação; separar inválidos, diagnósticos e viáveis com gap.
2. Ler base/auditoria antes das figuras; comparar custos e métricas, não só aparência.
3. Fronteira: 12 pontos, confirmações MILP e diferenças LP/MILP. Contornos são aproximados; veredito automático é preliminar, não substitui discussão de H2.
4. Comparar típico/P5–P15/P10: pico, economia FV, descarga/EFC BESS. NASA horária repetida em 15 min não prova o máximo mensal real.
5. Financeiro: dois semestres representativos, verão out–mar e inverno abr–set, calendário real. Segunda–sexta conta como dia de ponta; feriados só são excluídos se listados em `feriados_sem_ponta` num JSON de configuração (`--config arquivo.json`). Sem lista, a aproximação fica explícita.
6. Distinguir CAPEX* bruto, após regularizador e líquido FOM/reposição. Não somar o regularizador novamente ao pacote FOM. O&M FV não foi inventado: seus indicadores são brutos. FOM 2,5% é proxy internacional, não orçamento brasileiro.
7. Melhor contrato é o da **grade**, não ótimo contínuo. Com demanda faturável `max(D,P)`, contratos entre `P/1,05` e `P` podem empatar para o mesmo despacho.
8. Selecionar figuras para corpo/apêndice após avaliação. Nove disponíveis não significam nove obrigatórias. Conferir legibilidade no template.

Para retorno, compacte **somente a pasta da rodada**: já inclui código, parâmetros, tabelas e figuras. Não envie `.venv` nem replique o cache NASA. Nenhum e-mail foi enviado pelo pacote.
