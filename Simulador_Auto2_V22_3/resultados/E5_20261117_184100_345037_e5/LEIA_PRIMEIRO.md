# Leitura da rodada V22.3

EXECUÇÃO CONCLUÍDA

Modo: e5. Rodada de conferência não substitui a matriz E5.

Estados (incluem o mesmo caso em blocos distintos): {'CONVERGIDO': 277, 'HEURISTICA_VIAVEL': 77, 'HEURISTICA_COM_DEFICIT': 7, 'INVIAVEL_COMPROVADO': 3, 'VIAVEL_COM_GAP': 2}

1. Comece por tabelas/verificacao.csv e auditoria/auditoria.csv.

2. Em tabelas/cenarios.csv, filtre bloco=base: dez linhas nas duas modalidades.

3. VIAVEL_COM_GAP mantém a fatura do incumbente, sem afirmar ótimo. Consulte auditoria/solver.csv.

4. Limites fatura_min são do mínimo da fatura isolada. Em C3 não equivalem ao ótimo de fatura+regularizador.

5. financeiro_anual usa dois semestres representativos, não desempenho anual observado. Baixa geração é estresse.

6. fronteira usa CAPEX anual bruto e líquido; FOM já inclui reposição. Não somar ao regularizador.

7. As figuras têm células N/D e pontos inválidos omitidos explicitamente; não interpolar conclusões através de lacunas.

8. Verificar legibilidade no template antes de publicar. Fonte nominal 8 pt a 17,5 cm.

A execução do editor não foi realizada: esta rodada é evidência produzida pelo usuário.

H2_bruta: {'viaveis': 3, 'inviaveis': 9, 'indeterminadas': 0, 'veredito': 'corroborada nos pontos examinados', 'ressalva': 'grade aproximada; classes pontuais com gaps reportados; não certifica todo o contínuo entre pontos'}

H2_liquida: {'viaveis': 3, 'inviaveis': 9, 'indeterminadas': 0, 'veredito': 'corroborada nos pontos examinados', 'ressalva': 'grade aproximada; classes pontuais com gaps reportados; não certifica todo o contínuo entre pontos'}
