# V22.3 — alterações e limites

## Origem preservada

Derivada da V22.2 da etapa do e-mail ao professor, não da rodada posterior. Código e resultados anteriores permanecem intactos.

- V22.2: SHA256 `8b6f0be7e3421ece12501f3616ca335a0f8811e532651c3cee6395b89d26b254`.
- V22.2-r1 preservada: SHA256 `1a94047a0aa8819b331804a0384ff7fd7be3f6e7d12a23248ba9b818898a9d18`.
- Incorporada explicitamente a correção de gravação atômica com repetição em bloqueios do Windows.

## Causa → alteração → efeito esperado

| Causa/pedido | Alteração | Efeito a verificar |
|---|---|---|
| Arquivos excessivos | Tabelas consolidadas, registro único, dois CSVs de auditoria e PNGs | Rastreabilidade sem logs/arquivos por cenário |
| Sobrecarga com crédito zero | SoC final igual ao alvo; folga diagnóstica limitada | Mesma energia requerida; elimina consumo extra para absorver FV |
| Janela vazia | Igualdade como expressão PuLP mesmo sem variável de carga | Infactibilidade chega ao solver, sem igualdade booleana inválida |
| N=2 sem convergência | Melhor incumbente auditado e limites agregados entre tentativas | Viável com gap não é descartado como impossível |
| Limites da fatura | LB=max(0,LB_obj−R_max), UB=fatura do incumbente | Limites do mínimo da fatura isolada; não do ótimo fatura+regularizador |
| Warm start não transmitido | Envio explícito via highspy e registro de recebimento/fallback | Confirmar na rodagem; recebimento não prova uso/viabilidade |
| Tentativas demoradas | Teto acumulado por caso e progresso por bloco | Controle de duração preservando incumbentes |
| Base metodológica | Mantidos 54 veículos, 10 h, 300 kWh, 350 kW/carregador e parâmetros E4 | BESS não é obrigado a ciclar; 4 h continua sensibilidade |
| Transição de permanência | 4/6/8/10 h × crédito 0/0,5/1 | Examinar transição em 12 pontos |
| Fronteira | LP; MILP em extremos/base e vizinhos com mudança de classe bruta/líquida | Reportar divergências/lacunas, sem certificar domínio contínuo |
| Capacidade BESS | 500/750/1.000 kWh, 4 h, 500 kW fixos, créditos 0/1 | Comparar economia anual, CAPEX*, EFC, descarga e energia útil |
| Baixa geração | Mantida média P5–P15 e adicionado dia mais próximo de P10 | C2/C3 com mesma frota; data rastreável, sem alegar máximo mensal |
| Payback sazonal | Dois semestres ponderados por calendário real e 12 demandas | Estresse separado do financeiro anual |
| CAPEX líquido | Bruto, após regularizador e líquido FOM/reposição | Expor margem econômica sem duplicar degradação |
| Contrato | 400–800/50 kW; típico/P5–P15; C2/C3 | Ótimo da grade com demanda faturável=max(D,P) |
| Caso G | LP de arbitragem, perdas e extensão de ordem de grandeza | Não forçar igualdade entre arbitragem e economia total C2−C3 |
| Contas de William | Anuidade, CAPEX* e ultrapassagem na verificação | Conferências reproduzíveis com tolerância de arredondamento |
| Regularizador | 0,05 versus zero em pontos selecionados | Magnitude e proporção da economia; N/D se economia quase nula |
| Prazo E5 | Três sementes LP/MILP; dez configuráveis | Dispersão limitada à amostra examinada |
| Figuras | Até nove PNGs, nomes curtos, unidades e lacunas explícitas | Seleção para corpo/apêndice após avaliação |

## Ajustes posteriores no artigo

Igualdade de SoC é mais forte que o antigo mínimo e deve constar na formulação. Atualizar também grade 4×3, três sementes, confirmações MILP, calendário anual aproximado, P10 observado, contrato em grade e alternativas de CAPEX líquido. O artigo não foi editado nesta etapa.

Tarifas históricas E4 não foram atualizadas para orçamento comercial. Feriados, dois perfis anuais, NASA horária, crédito simplificado, ocupação livre e ausência de medições do operador continuam limitações. Não há uso mínimo, janela artificial nem benefício obrigatório do BESS.

## Fonte adicional do financeiro

ATB 2024b, **Commercial Battery Storage**, seção *Operation and Maintenance (O&M) Costs*: https://atb.nrel.gov/electricity/2024b/commercial_battery_storage . Consulta: 17/11/2026. FOM de 2,5% do capital inclui compensação da degradação para manter capacidade em 15 anos, com uso aproximadamente diário. A fração é aqui extrapolada como proxy de sensibilidade a 10 anos; não é cotação brasileira nem ajuste identificado para estes ciclos. Mantêm-se resultados bruto e após regularizador R$ 0,05/kWh. FOM/reposição e regularizador não são somados como custos físicos independentes.

## Revisão e pendências

Leitura estática do modelo, fluxo de arquivos e interfaces PuLP/highspy instaladas. Geração por texto/AST: simulador não foi importado nem executado. Análise sintática e busca de referências globais não substituem testes. Não foram produzidos resultados V22.3 nem renderizadas figuras novas. A rodagem do usuário deve confirmar funcionamento numérico, warm start, exportação e legibilidade. Não declarar validação antes disso.

Revisão científica continua após os resultados, especialmente com gaps, LP relaxado, confirmações incompletas ou economia comparável à incerteza do solver. O roteiro define a sequência e os critérios de retorno.

## Registro da entrega

SHA256 V22.3: `21bc1da390fe28dce18b34a17431519ad2c0e9e5bbaedd00cafd49a59c26d0fd`.

Revisao estatica: AST valido; referencias globais resolvidas; aridade de chamadas locais diretas conferida; hashes V22.2/r1 preservados; copia do cache identica. Nenhuma execucao/importacao do simulador.
