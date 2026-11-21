# -*- coding: utf-8 -*-
"""SIMEA 2027 — V22.3, revisão E5 após William.
Derivação textual da V22.2 (16/11), sem execução pelo editor.
Walz & Rudion: DOI 10.3390/en17225613. NASA POWER hourly 2019–2023.
Tarifas históricas e parâmetros físicos do E4; não orçamento de 2027.
Mudanças: alvo exato, incumbentes/limites, ano sazonal, fronteira 4×3,
P10 observado, capacidade BESS, três sementes, contrato 400–800/50 kW.
O&M/augmentation: ATB 2024b, Commercial Battery Storage, 2,5% CAPEX/ano:
https://atb.nrel.gov/electricity/2024b/commercial_battery_storage
É proxy internacional (referência de 15 anos), aqui extrapolada a 10 anos.
Não somar esse pacote à proxy de degradação de 0,05 R$/kWh: são ALTERNATIVAS.
Nunca se força ciclo/viabilidade do BESS. LP relaxado não é agenda física certificada.
"""
import os, sys, io, json, time, hashlib, platform, traceback, urllib.request
import argparse, calendar, copy, shutil
from pathlib import Path
from datetime import datetime, timezone, date
from importlib.metadata import version, PackageNotFoundError
try:
    import numpy as np
    import pandas as pd
    import pulp
    import highspy
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ModuleNotFoundError as erro:
    raise SystemExit(f'Dependência ausente: {erro.name}. Execute python -m pip install -r requirements.txt') from erro

VERSAO='22.3-E5'
CODIGO_ORIGEM_SHA256='8b6f0be7e3421ece12501f3616ca335a0f8811e532651c3cee6395b89d26b254'
DT=0.25
T_STEPS=96
IDX_PONTA=range(72,84)
_PASTA={'saida':'.'}
_ID_CHAMADA=[0]
REGISTRO_SOLVES=[]
REGISTRO_AUDITORIAS=[]
_CACHE_SOLAR_COMPLETO={}
_CACHE_FV={}
CONTEXTO_ATUAL=''
MEMORIA={}
REGISTRO={}
LINHAS=[]
AUDITORIAS=[]
COMPARACOES=[]
PERFIS=[]
VERIFICACOES=[]
FRONTEIRA=[]
FINANCEIRO=[]
DIMENSIONAMENTO=[]
FONTES_SOLARES={}


def _val(x):
    if hasattr(x, 'value'):
        v = x.value()
        if v is None:
            raise ValueError('Solver não forneceu valor para uma variável da solução.')
        return float(v)
    return float(x)


def _grade(vmin, vmax, passo, base):
    """Grade de varredura que SEMPRE contém o valor-base (evita ponto-base fora da matriz)."""
    return sorted(set(range(int(vmin), int(vmax) + 1, int(max(1, passo)))) | {int(base)})


def _finito(x):
    return x is not None and isinstance(x, (int, float, np.floating)) and np.isfinite(x)


def _extrair_gap(m, ret):
    """(gap relativo, bound dual) do solver, quando disponíveis; senão (nan, nan)."""
    gap, bound = float('nan'), float('nan')
    try:
        g = getattr(ret, 'gap_rel', None)                 # PuLP 4.x
        b = getattr(ret, 'best_bound', None)
        if g is not None:
            gap = float(g)
        if b is not None:
            bound = float(b)
    except Exception:
        pass
    if not np.isfinite(gap):
        try:
            info = m.solverModel.getInfo()                # PuLP 3.x + highspy
            gap = float(info.mip_gap)
            bound = float(info.mip_dual_bound)
        except Exception:
            pass
    if not np.isfinite(gap):
        gap = float('nan')
    return gap, bound


def baixar_nasa_power_horario(lat, lon, anos, tentar_dni_dif=True):
    """GHI, DNI e difusa horárias (W/m2, hora solar local LST), uma requisição por ano."""
    base = "https://power.larc.nasa.gov/api/temporal/hourly/point"
    pars = "ALLSKY_SFC_SW_DWN,ALLSKY_SFC_SW_DNI,ALLSKY_SFC_SW_DIFF" if tentar_dni_dif else "ALLSKY_SFC_SW_DWN"
    quadros = []
    for ano in anos:
        url = (f"{base}?parameters={pars}&community=RE&longitude={lon}&latitude={lat}"
               f"&start={ano}0101&end={ano}1231&format=CSV&time-standard=LST")
        with urllib.request.urlopen(url, timeout=180) as r:
            txt = r.read().decode('utf-8')
        linhas = txt.splitlines()
        i0 = next(i for i, l in enumerate(linhas) if l.startswith('YEAR'))
        quadros.append(pd.read_csv(io.StringIO("\n".join(linhas[i0:]))))
        time.sleep(1.0)
    df = pd.concat(quadros, ignore_index=True).replace(-999, np.nan)
    ren = {'ALLSKY_SFC_SW_DWN': 'GHI', 'ALLSKY_SFC_SW_DNI': 'DNI', 'ALLSKY_SFC_SW_DIFF': 'DIF'}
    df = df.rename(columns=ren)
    cols = ['YEAR', 'MO', 'DY', 'HR', 'GHI'] + (['DNI', 'DIF'] if tentar_dni_dif else [])
    return df[cols].reset_index(drop=True)


def irradiancia_plano_inclinado(df, lat, beta_graus, albedo=0.2):
    """G_i (W/m2) no plano inclinado voltado ao NORTE (hemisfério sul): direta (DNI x cos theta)
    + difusa isotrópica + refletida (Duffie & Beckman). DNI e DIF são obrigatórios; nenhuma transposição por média regional."""
    if 'DNI' not in df.columns or 'DIF' not in df.columns:
        raise ValueError('DNI/DIF obrigatórios na série científica V22.')
    datas = pd.to_datetime(dict(year=df['YEAR'], month=df['MO'], day=df['DY']))
    n = datas.dt.dayofyear.to_numpy(float)
    dec = np.radians(23.45 * np.sin(np.radians(360.0 * (284 + n) / 365.0)))
    w = np.radians(15.0 * (df['HR'].to_numpy(float) + 0.5 - 12.0))
    phi, b, g = np.radians(lat), np.radians(beta_graus), np.radians(180.0)   # g=180: voltado ao norte
    cos_t = (np.sin(dec) * np.sin(phi) * np.cos(b) - np.sin(dec) * np.cos(phi) * np.sin(b) * np.cos(g)
             + np.cos(dec) * np.cos(phi) * np.cos(b) * np.cos(w)
             + np.cos(dec) * np.sin(phi) * np.sin(b) * np.cos(g) * np.cos(w)
             + np.cos(dec) * np.sin(b) * np.sin(g) * np.sin(w))
    cos_z = np.sin(dec) * np.sin(phi) + np.cos(dec) * np.cos(phi) * np.cos(w)
    direta = np.where(cos_z > 0.02, df['DNI'].to_numpy(float) * np.clip(cos_t, 0, None), 0.0)
    difusa = df['DIF'].to_numpy(float) * (1 + np.cos(b)) / 2
    refletida = df['GHI'].to_numpy(float) * albedo * (1 - np.cos(b)) / 2
    return np.clip(direta + difusa + refletida, 0, None)


def gerar_pv(p, curva, pot_fv):
    pr = p['pr_fv'] if p.get('usar_pr_fv', True) else 1.0
    return np.repeat(curva, 4) * pot_fv * pr


def custo_kwh_mes_vetor(p, modalidade):
    """R$ mensais por kWh/dia importado em cada intervalo t, com o CALENDÁRIO da REN 1.000
    (Art. 2º, XXXVIII): a ponta só vale nos dias úteis. Mesma agenda diária todos os dias:
    dias_uteis x tarifa(t, com ponta) + (dias_mes - dias_uteis) x tarifa fora de ponta."""
    du, dm = p['dias_uteis_mes'], p['dias_mes']
    if modalidade == 'verde':
        pt, fp = p['tar_energia_ponta_verde'], p['tar_energia_fora_verde']
    else:
        pt, fp = p['tar_energia_ponta_azul'], p['tar_energia_fora_azul']
    return np.array([du * (pt if t in IDX_PONTA else fp) + (dm - du) * fp for t in range(T_STEPS)])


def tarifa_credito_injecao(p, modalidade):
    """Crédito (R$/kWh) da energia injetada: tarifa de energia fora de ponta menos a parcela NÃO
    compensada do Fio B (Lei 14.300/2022, Art. 27: 60% em 2026, 75% em 2027, 90% em 2028)."""
    fora = p['tar_energia_fora_verde'] if modalidade == 'verde' else p['tar_energia_fora_azul']
    return max(0.0, fora - p['pct_fio_b_nao_comp'] * p['fio_b_rs_kwh']) * p.get('multiplicador_credito', 1.0)


def custo_demanda_verde_fixo(pico, Dc, p):
    """Art. 301 da REN 1.000/2021 (contrato fixo, modalidade Verde). Retorna (demanda, ultrapassagem) em R$/mês.
    Ultrapassagem só se pico > (1+tol)*Dc; excedente = pico - Dc cobrado a mult x a tarifa.
    regra 'faturavel_max' (Eq. 3 do artigo): base = T*max(Dc, pico);  'excedente_2x': base = T*Dc.  regra adotada no E4"""
    tar, tol, mult = p['tar_dem_verde'], p['tolerancia_ultrapassagem'], p['mult_ultrapassagem']
    exc = max(0.0, pico - Dc) if pico > (1.0 + tol) * Dc else 0.0
    base = max(Dc, pico) * tar if p['regra_ultrapassagem'] == 'faturavel_max' else Dc * tar
    return base, mult * tar * exc


def fator_anuidade(i, n):
    return float(n) if i <= 0 else (1.0 - (1.0 + i) ** (-float(n))) / i


def gerar_chegadas_walz_rudion(f, seed):
    """DISTRIBUIÇÃO DE REFERÊNCIA: Walz e Rudion (2024, Energies, DOI 10.3390/en17225613), 'general cargo
    depot', Fig. 3b, blocos de 3 h (valores lidos do gráfico, não tabelados no artigo):
    0-3h 20% | 3-6h 20% | 15-18h 50% | 18-21h 10%."""
    rng = np.random.default_rng(seed)
    inicio = np.array([0, 3, 6, 9, 12, 15, 18, 21])
    prob = np.array([0.20, 0.20, 0.00, 0.00, 0.00, 0.50, 0.10, 0.00])
    return rng.choice(inicio, size=f, p=prob) + rng.uniform(0, 3, size=f)


def gerar_chegadas_noturno_hipotese(f, seed):
    """CENÁRIO DO AUTOR (hipótese, NÃO é dado medido): 60% 16h30-19h30 | 20% 19h30-22h30 | 20% 22h30-02h30."""
    rng = np.random.default_rng(seed)
    blocos = [(16.5, 3.0, 0.60), (19.5, 3.0, 0.20), (22.5, 4.0, 0.20)]
    idx = rng.choice(len(blocos), size=f, p=[b[2] for b in blocos])
    return np.array([(blocos[i][0] + rng.uniform(0.0, blocos[i][1])) % 24.0 for i in idx])


def gerar_frota(f, p, seed=42):
    """Perfil individual de cada caminhão: chegada, janela, SoC inicial, capacidade, energia/turno.
    Janela arredondada PARA DENTRO: chegada em ceil e partida em floor (múltiplos de 15 min).
    Frotas ANINHADAS: a frota de f veículos é o prefixo da frota de referência (mesma semente), de modo
    que cenários e tamanhos de frota comparados usam as MESMAS realizações."""
    noturno = p['perfil_chegada'] == 'noturno_hipotese'
    fr = max(int(f), int(p.get('frota_max', f)), int(p.get('frota_base') or f))
    cheg_h = gerar_chegadas_noturno_hipotese(fr, seed) if noturno else gerar_chegadas_walz_rudion(fr, seed)
    rng = np.random.default_rng(seed + 10_000)
    soc_i = np.clip(rng.normal(p['soc_i_mean'], p['soc_i_std'], fr), p['soc_i_min'], p['soc_i_max'])
    perm_sorteada = np.clip(rng.normal(p['perm_media_h'], p['perm_desvio_h'], fr), p['perm_min_h'], p['perm_max_h'])
    if noturno and p.get('limitar_partida_manha', True):
        limite = p['hora_limite_partida_manha'] + rng.uniform(-0.5, 0.5, fr)
        ate_manha = (limite - cheg_h) % 24.0
        # a permanência efetiva pode ficar abaixo do mínimo nominal (partida matinal manda)
        perm_h = np.maximum(np.minimum(perm_sorteada, ate_manha), DT)
    else:
        perm_h = perm_sorteada
    cheg_h, perm_h, soc_i = cheg_h[:f], perm_h[:f], soc_i[:f]
    a_abs = np.ceil(cheg_h / DT - 1e-9).astype(int)
    fim_abs = np.floor((cheg_h + perm_h) / DT + 1e-9).astype(int)
    L = np.clip(fim_abs - a_abs, 0, T_STEPS)  # janela vazia gera diagnóstico, nunca intervalo inventado
    a = a_abs % T_STEPS
    cap = np.full(f, float(p['cap_bat']))
    return {'f': f, 'a': a, 'L': L, 'cap': cap, 'soc_i': soc_i, 'e_nec': cap * (p['soc_f'] - soc_i),
            'chegada_h': cheg_h, 'partida_h': ((a_abs + L) * DT) % 24.0, 'perm_h': L * DT}


def pmax_carregador(soc, p):
    """Potência máxima admissível (kW, lado CA) em função do SoC no início do intervalo: P_c até
    soc_cv e redução linear até zero em SoC = 1 (aproximação CC/CV). Usada em C0, no diagnóstico de
    janelas, na auditoria e (de forma equivalente) nas restrições do MILP."""
    s = np.asarray(soc, dtype=float)
    pm = p['pot_carregador']
    if p['usar_ccv']:
        return np.minimum(pm, pm * (1.0 - s) / (1.0 - p['soc_cv']))
    return pm * np.ones_like(s)


def simular_nao_coordenada(frota, p, n_carregadores=None, bloquear_ponta=False):
    """Cada caminhão conecta ao chegar se houver carregador livre (FIFO), carrega à potência máxima
    (P_c e curva CC/CV, a MESMA função pmax_carregador do EMS) até o SoC alvo ou a partida; libera o
    carregador ao terminar. Regime periódico: simulam-se 8 dias e usa-se o último, verificando repetição periódica. Devolve também o
    cronograma POR VEÍCULO dobrado no dia (P_veic), usado no teste de coerência física e como partida
    a quente do C1. bloquear_ponta=True (C0b): regra simples, o veículo espera plugado até o fim da ponta."""
    N = n_carregadores or p['n_carregadores']
    f, n_dias = frota['f'], 8
    eta = p['eta_carga']
    cheg = np.concatenate([frota['a'] + T_STEPS * d for d in range(n_dias)])
    part = cheg + np.tile(frota['L'], n_dias)
    capv = np.tile(frota['cap'], n_dias)
    e = np.tile(frota['soc_i'] * frota['cap'], n_dias).astype(float)
    alvo = capv * p['soc_f']
    nv = len(cheg)
    ordem = np.argsort(cheg, kind='stable')
    inicio = np.full(nv, np.nan)
    carga = np.zeros(n_dias * T_STEPS)
    Pinst = np.zeros((nv, n_dias * T_STEPS))
    esperando, plugado, prox = [], [], 0
    for s in range(n_dias * T_STEPS):
        while prox < nv and cheg[ordem[prox]] <= s:
            esperando.append(ordem[prox])
            prox += 1
        esperando = [i for i in esperando if part[i] > s]
        plugado = [i for i in plugado if part[i] > s]
        while len(plugado) < N and esperando:
            i = esperando.pop(0)
            plugado.append(i)
            inicio[i] = s
        restantes = []
        for i in plugado:
            soc = e[i] / capv[i]
            pmax = float(pmax_carregador(soc, p))
            pi = max(0.0, min(pmax, (alvo[i] - e[i]) / (eta * DT)))
            if bloquear_ponta and (s % T_STEPS) in IDX_PONTA:
                pi = 0.0
            carga[s] += pi
            Pinst[i, s] = pi
            e[i] += eta * DT * pi
            if e[i] < alvo[i] - 1e-6:
                restantes.append(i)
        plugado = restantes
    meio = np.arange((n_dias - 2) * f, (n_dias - 1) * f)
    inicio_m = np.where(np.isnan(inicio[meio]), part[meio], inicio[meio])
    espera_h = (inicio_m - cheg[meio]) * DT
    deficit = np.maximum(0.0, alvo[meio] - e[meio])
    deficit = np.where(deficit < 1e-3, 0.0, deficit)
    P_veic = np.zeros((f, T_STEPS))
    for d in range(n_dias):
        P_veic += Pinst[d * f:(d + 1) * f, (n_dias - 1) * T_STEPS:n_dias * T_STEPS]
    return {'carga': carga[(n_dias - 1) * T_STEPS:n_dias * T_STEPS].copy(), 'P_veic': P_veic,
            'residuo_periodico_kw': float(np.max(np.abs(carga[-T_STEPS:] - carga[-2*T_STEPS:-T_STEPS]))),
            'espera_media_h': float(np.mean(espera_h)), 'espera_maxima_h': float(np.max(espera_h)),
            'pct_sem_espera': float(np.mean(espera_h < 1e-6) * 100),
            'energia_nao_atendida': float(np.sum(deficit)), 'n_nao_atendidos': int(np.sum(deficit > 0))}


def _c_deg(p):
    """Coeficiente de desgaste usado NO OBJETIVO (R$/kWh descarregado). Padrão: regularizador pequeno."""
    v = p['custo_degrad_otim']
    return resolver_lcos(p) if isinstance(v, str) else float(v)


def otimizar_ems(frota, p, cap_bess, pot_fv, curva, modalidade='verde', modo=None,
                 n_carregadores=None, dem_contratada=None, rotulo='', warm=None,
                 permitir_deficit=False, tempo_limite=None, tentativa=1, id_chamada=0):
    """
    Conjuntos: v (veículos), t = 0..95 (15 min, cíclico), k = 0..L_v-1 (janela do veículo v).
    Variáveis: p[v,k] (kW), soc_v[v,k] in [0,1], u[v,k] in {0,1} (carregador ocupado; MILP), ini[v,k]
    (só ocupacao='continua'), P_imp[t], P_exp[t], P_ch[t], P_dis[t] (BESS, variáveis separadas, >= 0),
    soc_b[t], P_pico (Verde) ou P_pico_p/P_pico_f (Azul) por epígrafe, C_energia_liq >= 0; no contrato
    fixo: z (binária de ultrapassagem), U_exc, D_fat. Só no DIAGNÓSTICO: folga s[v] (kWh) por veículo.
    Restrições: balanço (Eq. 1/5), SoC dinâmico por veículo com rendimento, janela, SoC ALVO DURO
    (soc_v[L] == soc_f), potência do carregador (+ CC/CV), simultaneidade (sum u <= N), limites do BESS
    e fechamento periódico do ciclo, epígrafe do pico. SEM janela horária imposta ao BESS.
    Objetivo: C_demanda + C_energia_liq + regularizador de desgaste (+ penalidade de folga no diagnóstico).
    """
    t0 = time.time()
    modo = (modo or p['modo_simultaneidade']).lower()
    N = n_carregadores or p['n_carregadores']
    D = dem_contratada if dem_contratada is not None else (
        p['dem_contratada'] if p.get('contrato_modo', 'otimo') == 'fixo' else None)
    if modalidade != 'verde':
        D = None     # contrato fixo (Art. 301) só implementado para a modalidade Verde
    continua = p.get('ocupacao', 'livre') == 'continua'
    f, a, L, cap = frota['f'], frota['a'], frota['L'], frota['cap']
    P_c, eta, eta_b = p['pot_carregador'], p['eta_carga'], p['eta_bess']
    DM = p['dias_mes']
    pv_gen = gerar_pv(p, curva, pot_fv)
    w_en = custo_kwh_mes_vetor(p, modalidade)
    tar_dem_v = p['tar_dem_verde']
    tem_bess = cap_bess > 0
    c_deg = _c_deg(p)
    tl = tempo_limite if tempo_limite else p['tempo_limite_s']

    presenca = np.zeros(T_STEPS, dtype=int)
    for v in range(f):
        for k in range(L[v]):
            presenca[(a[v] + k) % T_STEPS] += 1
    simult_ativa = presenca.max() > N
    usar_bin = simult_ativa and modo == 'milp'

    m = pulp.LpProblem("EMS_por_veiculo", pulp.LpMinimize)
    cont_bin, cont_var = [0], [0]

    def var(nome, lo=None, hi=None, binaria=False):
        cont_var[0] += 1
        if binaria:
            cont_bin[0] += 1
        cat = 'Binary' if binaria else 'Continuous'
        if hasattr(m, 'add_variable'):
            return m.add_variable(nome, lowBound=lo, upBound=hi, cat=cat)
        return pulp.LpVariable(nome, lowBound=lo, upBound=hi, cat=cat)

    pv_var, u_var, socv = {}, {}, {}
    carga_t = [[] for _ in range(T_STEPS)]
    u_t = [[] for _ in range(T_STEPS)]
    folga = {}
    for v in range(f):
        soc0 = float(frota['soc_i'][v])
        soc_prev = soc0
        kcoef = P_c / (1.0 - p['soc_cv'])
        u_prev, inicios = None, []
        for k in range(L[v]):
            pk = var(f"p_{v}_{k}", 0, P_c)
            pv_var[(v, k)] = pk
            t = (a[v] + k) % T_STEPS
            carga_t[t].append(pk)
            if usar_bin:
                uk = var(f"u_{v}_{k}", binaria=True)
                u_var[(v, k)] = uk
                m += pk <= P_c * uk
                u_t[t].append(uk)
                if continua:
                    ini = var(f"ini_{v}_{k}", 0, 1)
                    m += (ini >= uk) if u_prev is None else (ini >= uk - u_prev)
                    inicios.append(ini)
                    u_prev = uk
            if p['usar_ccv'] and (soc0 + (eta * DT * P_c * k) / cap[v] > p['soc_cv']):
                m += pk + kcoef * soc_prev <= kcoef
            nxt = var(f"soc_{v}_{k + 1}", 0.0, 1.0)
            m += nxt == soc_prev + (eta * DT * pk) / float(cap[v])
            socv[(v, k + 1)] = nxt
            soc_prev = nxt
        if usar_bin and continua:
            m += pulp.lpSum(inicios) <= 1
        if permitir_deficit:                       # DIAGNÓSTICO de infactibilidade (folga em kWh)
            s_v = var(f"s_{v}", 0, max(0., float(cap[v])*(p['soc_f']-soc0)))
            folga[v] = s_v
            m += soc_prev + s_v / float(cap[v]) == float(p['soc_f'])
        else:                                      # SoC ALVO = RESTRIÇÃO DURA
            m += pulp.lpSum([soc_prev]) == float(p['soc_f'])

    if simult_ativa:
        for t in range(T_STEPS):
            if presenca[t] > N:
                m += (pulp.lpSum(u_t[t]) <= N) if usar_bin else (pulp.lpSum(carga_t[t]) <= N * P_c)

    p_imp = [var(f"P_imp_{t}", 0) for t in range(T_STEPS)]
    p_exp = [var(f"P_exp_{t}", 0) for t in range(T_STEPS)]
    zb_var = []
    if tem_bess:
        p_ch = [var(f"P_ch_{t}", 0, p['pot_bess']) for t in range(T_STEPS)]
        p_dis = [var(f"P_dis_{t}", 0, p['pot_bess']) for t in range(T_STEPS)]
        soc_b = [var(f"SoC_b_{t}", p['soc_min_bess'], p['soc_max_bess']) for t in range(T_STEPS)]
        if p.get('bess_exclusao', 'posteriori') == 'binaria' and modo == 'milp':
            for t in range(T_STEPS):
                zb = var(f"zb_{t}", binaria=True)
                zb_var.append(zb)
                m += p_ch[t] <= p['pot_bess'] * zb
                m += p_dis[t] <= p['pot_bess'] * (1 - zb)
    for t in range(T_STEPS):
        ch = p_ch[t] if tem_bess else 0
        dis = p_dis[t] if tem_bess else 0
        m += pv_gen[t] + p_imp[t] + dis == pulp.lpSum(carga_t[t]) + p_exp[t] + ch          # Eq. (1)/(5)
        if tem_bess:   # SoC cíclico: soc_b[t] ligado a soc_b[t-1], com t-1 = T-1 para t = 0 (ciclo periódico fechado)
            m += soc_b[t] == soc_b[(t - 1) % T_STEPS] + (eta_b * p_ch[t] * DT - (p_dis[t] * DT) / eta_b) / float(cap_bess)

    big_m = float(p['big_m_kw'])
    pico_vars = {}
    z_ult = exc_ult = d_fat = None
    if modalidade == 'verde':
        p_pico = var("P_pico", 0)
        pico_vars['verde'] = p_pico
        for t in range(T_STEPS):
            m += p_pico >= p_imp[t]
        if D is not None:
            z_ult = var("z_ult", binaria=True)
            exc_ult = var("U_exc", 0)
            m += p_pico <= (1.0 + p['tolerancia_ultrapassagem']) * D + big_m * z_ult
            m += exc_ult >= p_pico - D - big_m * (1 - z_ult)
            if p['regra_ultrapassagem'] == 'faturavel_max':
                d_fat = var("D_fat", D)
                m += d_fat >= p_pico
                custo_dem = tar_dem_v * d_fat + p['mult_ultrapassagem'] * tar_dem_v * exc_ult
            else:
                custo_dem = tar_dem_v * D + p['mult_ultrapassagem'] * tar_dem_v * exc_ult
        else:
            custo_dem = tar_dem_v * p_pico
    else:
        pp = var("P_pico_p", 0)
        pf = var("P_pico_f", 0)
        pico_vars['ponta'], pico_vars['fora'] = pp, pf
        for t in range(T_STEPS):
            m += (pp if t in IDX_PONTA else pf) >= p_imp[t]
        custo_dem = p['tar_dem_ponta_azul'] * pp + p['tar_dem_fora_azul'] * pf
    tar_cred = tarifa_credito_injecao(p, modalidade)

    custo_energia_bruta = pulp.lpSum([p_imp[t] * float(w_en[t]) * DT for t in range(T_STEPS)])
    credito = pulp.lpSum([p_exp[t] * tar_cred * DT * DM for t in range(T_STEPS)])
    c_liq = var("C_energia_liq", 0)       # crédito nunca torna a fatura negativa
    m += c_liq >= custo_energia_bruta - credito
    desgaste = pulp.lpSum([p_dis[t] * DT * c_deg * DM for t in range(T_STEPS)]) if tem_bess else 0
    penal = p['custo_deficit_diag'] * pulp.lpSum(folga.values()) if permitir_deficit else 0
    m += custo_dem + c_liq + desgaste + penal

    # ---------- partida a quente (melhor esforço; nunca interrompe a otimização) ----------
    warm_ok = False
    usar_warm = (warm is not None and p.get('warm_start', True) and usar_bin and not continua
                 and not permitir_deficit and warm.get('P_veic') is not None)
    if usar_warm:
        try:
            P0 = warm['P_veic']
            vals = {}
            load = np.zeros(T_STEPS)
            for v in range(f):
                s_acc = float(frota['soc_i'][v])
                for k in range(L[v]):
                    t = (a[v] + k) % T_STEPS
                    pk = float(min(max(P0[v, t], 0.0), P_c))
                    load[t] += pk
                    vals[pv_var[(v, k)]] = pk
                    if (v, k) in u_var:
                        vals[u_var[(v, k)]] = 1.0 if pk > 1e-6 else 0.0
                    s_acc = min(1.0, s_acc + eta * DT * pk / float(cap[v]))
                    vals[socv[(v, k + 1)]] = s_acc
            ch0 = np.array(warm.get('bess_in', np.zeros(T_STEPS)), dtype=float) if tem_bess else np.zeros(T_STEPS)
            dis0 = np.array(warm.get('bess_out', np.zeros(T_STEPS)), dtype=float) if tem_bess else np.zeros(T_STEPS)
            net = load + ch0 - dis0 - pv_gen
            imp0, exp0 = np.maximum(net, 0.0), np.maximum(-net, 0.0)
            for t in range(T_STEPS):
                vals[p_imp[t]], vals[p_exp[t]] = float(imp0[t]), float(exp0[t])
            if tem_bess:
                if 'bess_soc' in warm:
                    socb0 = np.clip(np.array(warm['bess_soc'], dtype=float) / 100.0, p['soc_min_bess'], p['soc_max_bess'])
                else:
                    socb0 = np.full(T_STEPS, p['soc_min_bess'])
                for t in range(T_STEPS):
                    vals[p_ch[t]], vals[p_dis[t]], vals[soc_b[t]] = float(ch0[t]), float(dis0[t]), float(socb0[t])
                for t, zb in enumerate(zb_var):
                    vals[zb] = 1.0 if ch0[t] > 1e-6 else 0.0
            if modalidade == 'verde':
                pico0 = float(imp0.max())
                vals[pico_vars['verde']] = pico0
                if D is not None:
                    zz = 1.0 if pico0 > (1.0 + p['tolerancia_ultrapassagem']) * D else 0.0
                    vals[z_ult] = zz
                    vals[exc_ult] = max(0.0, pico0 - D) if zz > 0 else 0.0
                    if d_fat is not None:
                        vals[d_fat] = max(D, pico0)
            else:
                vals[pico_vars['ponta']] = float(imp0[list(IDX_PONTA)].max())
                vals[pico_vars['fora']] = float(max(imp0[:72].max(), imp0[84:].max()))
            bruta0 = float(np.sum(imp0 * w_en * DT))
            cred0 = float(np.sum(exp0 * tar_cred * DT * DM))
            vals[c_liq] = max(0.0, bruta0 - cred0)
            if len(vals) == cont_var[0]:
                for vr, val in vals.items():
                    vr.setInitialValue(float(val), check=False)
                warm_ok = True
        except Exception:
            warm_ok = False

    # ---------- resolver ----------
    solver, nome_solver = criar_solver(p, tl, warm_ok)
    try:
        ret = m.solve(solver)
    except Exception:
        if not warm_ok:
            raise
        warm_ok = False
        solver, nome_solver = criar_solver(p, tl, False)
        ret = m.solve(solver)
    dt_solve = time.time() - t0
    status, tem_sol = _status_modelo(m, ret)
    gap, bound = _extrair_gap(m, ret)
    tem_bin = cont_bin[0] > 0
    # PuLP 3: status Optimal sozinho pode ser incumbente em limite de tempo.
    # Fonte: https://coin-or.github.io/pulp/technical/constants.html
    sol_status = getattr(m, 'sol_status', None)
    otimo_reportado = (status == 'Optimal' and
                       (sol_status == getattr(pulp, 'LpSolutionOptimal', 1)))
    convergido = bool(tem_sol and (
        (gap <= p['gap_rel'] + 1e-9) if tem_bin and np.isfinite(gap) else otimo_reportado))
    if not tem_bin and otimo_reportado:
        gap, bound = 0.0, float(pulp.value(m.objective))
    infeasivel = (status == 'Infeasible')
    if not tem_sol and not infeasivel:
        print(f"   [AVISO] {rotulo} sem solução utilizável (status: {status}, {dt_solve:.0f} s).")
    elif tem_sol and not convergido:
        print(f"   [AVISO] {rotulo} solução viável sem fechar o gap (status {status}, gap "
              f"{gap * 100 if np.isfinite(gap) else float('nan'):.2f}%, {dt_solve:.0f} s).")

    try:
        n_var = len(m.variables())
    except Exception:
        n_var = cont_var[0]
    obj = float('nan')
    if tem_sol:
        try:
            obj = float(pulp.value(m.objective))
        except Exception:
            obj = float('nan')
    REGISTRO_SOLVES.append({'id': id_chamada, 'tentativa': tentativa, 'rotulo': rotulo,
                            'modo': ('MILP' if tem_bin else 'LP'), 'ocupacao_relaxada': simult_ativa and not usar_bin,
                            'solver': nome_solver, 'warm_preparado': warm_ok,
                            'warm_solicitado': bool(getattr(solver, 'inicio_enviado', False)),
                            'warm_aceito': getattr(solver, 'inicio_status', 'não solicitado'), 'tempo_s': dt_solve, 'limite_tempo_s': tl, 'f': f,
                            'n_bin': cont_bin[0], 'n_var': n_var, 'status': status,
                            'gap_pct': gap * 100 if np.isfinite(gap) else float('nan'),
                            'convergido': convergido, 'infeasivel': infeasivel, 'obj': obj,
                            'diagnostico': permitir_deficit, 'bound': bound, 'sol_status': sol_status})
    base = {'status': status, 'tem_solucao': bool(tem_sol), 'infeasivel': infeasivel, 'convergido': convergido,
            'gap': gap, 'bound': bound, 'tempo_s': dt_solve, 'n_bin': cont_bin[0], 'warm_ok': warm_ok,
            'tentativa': tentativa, 'obj': obj, 'diagnostico': permitir_deficit,
            'solver': nome_solver, 'sol_status': sol_status, 'otimo_reportado': otimo_reportado,
            'erro_objetivo_estimado_rs': max(0.0, obj-bound) if np.isfinite(obj) and np.isfinite(bound) else np.nan,
            'ocupacao_relaxada': simult_ativa and not usar_bin, 'N_aud': N,
            'warm_preparado': warm_ok, 'warm_solicitado': bool(getattr(solver, 'inicio_enviado', False)), 'warm_aceito': getattr(solver, 'inicio_status', 'não solicitado')}
    if not tem_sol:
        return base

    P_veic = np.zeros((f, T_STEPS))
    for (v, k), vr in pv_var.items():
        P_veic[v, (a[v] + k) % T_STEPS] += _val(vr)
    imp = np.array([_val(p_imp[t]) for t in range(T_STEPS)])
    exp_ = np.array([_val(p_exp[t]) for t in range(T_STEPS)])
    ch_ = np.array([_val(p_ch[t]) for t in range(T_STEPS)]) if tem_bess else np.zeros(T_STEPS)
    dis_ = np.array([_val(p_dis[t]) for t in range(T_STEPS)]) if tem_bess else np.zeros(T_STEPS)
    soc_bess = np.array([_val(soc_b[t]) for t in range(T_STEPS)]) * 100.0 if tem_bess else np.zeros(T_STEPS)
    recarga = P_veic.sum(axis=0)
    deficit_folga = sum(_val(s) for s in folga.values())
    residuo = np.max(np.abs(pv_gen + imp + dis_ - recarga - exp_ - ch_))
    soc_final = frota['soc_i'] + (eta * DT * P_veic.sum(axis=1)) / cap
    deficit = float(np.maximum(0.0, (p['soc_f'] - soc_final) * cap).sum())
    n_def = int(np.sum((p['soc_f'] - soc_final) * cap > 1e-3))
    simult = int(np.max((P_veic > 1e-6).sum(axis=0))) if f > 0 else 0
    ocup = np.zeros(T_STEPS)
    for (v, k), uvar in u_var.items():
        ocup[(a[v] + k) % T_STEPS] += _val(uvar)
    max_ocup = int(round(ocup.max())) if usar_bin else simult
    base.update({'grid_import': imp, 'grid_export': exp_, 'bess_in': ch_, 'bess_out': dis_, 'bess_soc': soc_bess,
                 'recarga': recarga, 'P_veic': P_veic, 'pv_gen': pv_gen, 'energia_nao_atendida': deficit,
                 'n_nao_atendidos': n_def, 'residuo_balanco_kw': residuo, 'max_simultaneos': simult,
                 'max_ocupacao': max_ocup, 'deficit_folga_kwh': deficit_folga, 'presenca': presenca, 'soc_final': soc_final,
                 'max_min_bess': float(np.max(np.minimum(ch_, dis_))),
                 'max_min_impexp': float(np.max(np.minimum(imp, exp_)))})
    return base


def calcular_fatura(imp, exp_, p, modalidade, dem_contratada):
    """Fatura mensal (dias_mes dias; ponta só nos dias úteis). OPEX = demanda + energia líquida +
    ultrapassagem, SEM qualquer custo fixo. Contrato 'otimo' (dem_contratada=None): demanda = pico.
    Contrato fixo (só Verde): regra do Art. 301 via custo_demanda_verde_fixo (a mesma do otimizador)."""
    w = custo_kwh_mes_vetor(p, modalidade)
    energia_bruta = float(np.sum(imp * w * DT))
    credito = float(np.sum(exp_ * tarifa_credito_injecao(p, modalidade) * DT) * p['dias_mes'])
    energia_liq = max(0.0, energia_bruta - credito)
    pico = float(np.max(imp))
    if modalidade == 'verde':
        if dem_contratada is None:
            dem, ult = pico * p['tar_dem_verde'], 0.0
        else:
            dem, ult = custo_demanda_verde_fixo(pico, dem_contratada, p)
    else:
        pico_p = float(np.max(imp[list(IDX_PONTA)]))
        pico_f = float(max(np.max(imp[:72]), np.max(imp[84:])))
        dem, ult = pico_p * p['tar_dem_ponta_azul'] + pico_f * p['tar_dem_fora_azul'], 0.0
    return {'demanda': dem, 'ultrapassagem': ult, 'energia': energia_liq,
            'energia_bruta': energia_bruta, 'credito': credito, 'total': dem + ult + energia_liq}


def avaliar(imp, exp_, pv_gen, bess_out, cap_bess, p, modalidade, dem_contratada, extras=None):
    pico = float(np.max(imp))
    fat = calcular_fatura(imp, exp_, p, modalidade, dem_contratada)
    total_fv = float(np.sum(pv_gen * DT))
    inj = float(np.sum(exp_ * DT))
    autoc = (total_fv - inj) / total_fv * 100.0 if total_fv > 0 else float('nan')
    e_dis_dia = float(np.sum(bess_out * DT))
    r = {'pico': pico, 'fator_carga': float(np.mean(imp)) / pico if pico > 0 else float('nan'),
         'custo_demanda': fat['demanda'], 'custo_ultrapassagem': fat['ultrapassagem'],
         'custo_energia': fat['energia'], 'energia_bruta_rs': fat['energia_bruta'],
         'credito_potencial_rs': fat['credito'], 'credito_utilizado_rs': min(fat['credito'], fat['energia_bruta']),
         'regularizador_rs_mes': e_dis_dia * p['dias_mes'] * _c_deg(p),
         'importacao_kwh_mes': float(np.sum(imp) * DT * p['dias_mes']),
         'ponta_kwh_mes': float(np.sum(imp[list(IDX_PONTA)]) * DT * p['dias_uteis_mes']),
         'opex': fat['total'], 'autoconsumo': autoc,
         'efc': e_dis_dia * p['dias_mes'] / cap_bess if cap_bess > 0 else 0.0,     # ciclos equivalentes por MÊS
         'efc_dia': e_dis_dia / cap_bess if cap_bess > 0 else 0.0,
         'degradacao_lcos': e_dis_dia * p['dias_mes'] * resolver_lcos(p),           # R$/mês, indicador ex-post
         'injecao_kwh_dia': inj, 'ponta_kwh': float(np.sum(imp[list(IDX_PONTA)]) * DT), 'imp': imp, 'exp': exp_}
    if extras:
        r.update(extras)
    return r


def obter_input(texto_prompt, valor_padrao, tipo=float, minusculo=False):
    entrada = input(f"{texto_prompt} [{valor_padrao}]: ").strip()
    if not entrada:
        return valor_padrao
    try:
        if tipo == str:
            return entrada.lower() if minusculo else entrada
        return tipo(entrada)
    except ValueError:
        print(f"Erro de digitação. Mantendo valor padrão: {valor_padrao}")
        return valor_padrao


def parametros_padrao():
    p = {
        # --- Localidade e série solar ---
        'estacao': 'Inverno',
        'lat_graus': -23.45, 'lon_graus': -46.53,       # Guarulhos-SP
        'inclinacao_fv_graus': 23.0,                    # ≈ latitude, módulos voltados ao norte
        'pr_fv': 0.80, 'usar_pr_fv': True,              # Atlas INPE (2017), p. 57: taxa de desempenho de 80%
        'arquivo_serie_fv': 'dados_nasa_v22',
        'baixar_serie_fv': True, 'anos_serie_fv': [2019, 2020, 2021, 2022, 2023],
        'area_fv_m2_por_kwp': 6.0,                      # [PREMISSA] só para o aviso de área

        # --- Frota e carregadores ---
        # 537 cam./ha/dia = FLUXO diário do centro logístico de W&R (Tab. 2); 0,10 ha é a área de REFERÊNCIA
        # que dimensiona o fluxo (54 veículos/turno), NÃO o terreno onde o FV é instalado.
        'area_ref_fluxo_ha': 0.10, 'caminhoes_ha_dia_cargo_geral': 537,
        'frota_base': 54,                             # None => round(area_ref x 537) = 54
        'frota_min': 20, 'frota_max': 100, 'frota_passo': 20,
        'n_carregadores': 10, 'carregadores_lista': [1, 2, 3, 4, 6, 8, 10],
        'pot_carregador': 350.0, 'cap_bat': 300.0,      # escolhas E4: faixa Fig.4a e cenário Fig.11, Walz/Rudion DOI 10.3390/en17225613; não MCS
        'soc_i_mean': 0.20, 'soc_i_std': 0.05, 'soc_i_min': 0.05, 'soc_i_max': 0.40, 'soc_f': 0.90,
        'eta_carga': 0.95,
        'usar_ccv': True, 'soc_cv': 0.80,
        'perfil_chegada': 'walz_rudion',                # CASO-BASE: 'walz_rudion' (referência da literatura) | 'noturno_hipotese' (autor)
        'limitar_partida_manha': True,                  # só afeta o perfil noturno_hipotese
        'hora_limite_partida_manha': 7.5,
        'perm_media_h': 10.0, 'perm_desvio_h': 2.0, 'perm_min_h': 6.0, 'perm_max_h': 14.0,   # permanência: HIPÓTESE
        'ocupacao': 'livre',                            # 'livre': pausa libera o carregador | 'continua': 1 bloco contíguo

        # --- Microrrede ---
        'pot_fv': 1200.0,
        'cap_bess': 1000.0, 'pot_bess': 500.0, 'soc_min_bess': 0.10, 'soc_max_bess': 1.00, 'eta_bess': 0.95,
        'ciclos_vida_util_bess': None,
        'bess_exclusao': 'posteriori',                  # 'posteriori' (verifica) | 'binaria' (impõe c.d = 0; mais lento)

        # --- Calendário e tarifas EDP SP A4, REH ANEEL 3.541/2025 — valores históricos aprovados no E4 ---
        'dias_mes': 30, 'dias_uteis_mes': 22,           # ponta só nos dias úteis (REN 1.000, Art. 2º, XXXVIII)
        'contrato_modo': 'otimo', 'dem_contratada': 2000.0, 'tolerancia_ultrapassagem': 0.05,
        'mult_ultrapassagem': 2.0,
        'regra_ultrapassagem': 'faturavel_max',         # 'faturavel_max' (Eq. 3) | 'excedente_2x' regra do E4, arts. 294 e 301
        'big_m_kw': 20000.0,
        'tar_dem_verde': 12.02, 'tar_dem_ponta_azul': 48.65, 'tar_dem_fora_azul': 12.02,
        'tar_energia_ponta_verde': 1.81466, 'tar_energia_fora_verde': 0.45128,
        'tar_energia_ponta_azul': 0.63108, 'tar_energia_fora_azul': 0.45128,   # FP Verde e Azul iguais na base documentada do E4
        'pct_fio_b_nao_comp': 0.75,                     # Lei 14.300, Art. 27: 75% em 2027
        'fio_b_rs_kwh': 0.10,                           # hipótese explícita do E4, não tarifa homologada
        'multiplicador_credito': 1.0,                   # 1,0 = crédito-base simplificado (caso-base); 0,0 = injeção sem crédito

        # --- CAPEX e break-even (CAPEX e OPEX comparados à parte; nada de CAPEX dentro do OPEX) ---
        'custo_pv_kwp': 2700.0,
        'custo_bess_kwh': 1500.0,                       # CAPEX de REFERÊNCIA do BESS (define a fronteira)
        'capex_bess_faixa': [1100.0, 1500.0, 3000.0, 4000.0],   # R$/kWh: recorte autoral da faixa explorada pela EPE, conforme E4
        'capex_fv_faixa': [2700.0, 3630.0],             # R$/kWp
        'taxa_desconto': 0.10, 'vida_bess_anos': 10, 'vida_fv_anos': 25,
        'lcos_bess': None,                            # desativado: sem hipótese de ciclos de vida útil não usada no E4
        'custo_degrad_otim': 0.05,                      # R$/kWh descarregado NO OBJETIVO (regularizador pequeno)
        'custo_deficit_diag': 1000.0,                   # R$/kWh: folga SÓ no diagnóstico de infactibilidade (>> custo da energia)

        # --- Solver: protocolo de convergência ---
        'solver': 'HiGHS', 'tempo_limite_s': 60, 'tempo_max_s': 600, 'escalonar_tempo': True,
        'gap_rel': 0.002, 'warm_start': True, 'solver_msg': False,
        'modo_simultaneidade': 'milp', 'modo_varredura': 'lp',   # grade da fronteira / Monte Carlo / contrato em LP relaxado
        'n_realizacoes_monte_carlo': 3, 'rodar_verificacao': True, 'modo_rapido': False,

        # --- Varreduras ---
        'perm_sweep_h': [4.0, 6.0, 8.0, 10.0],    # eixo permanência (Tabela 4 e fronteira)
        'credito_sweep': [0.0, 0.5, 1.0],   # eixo crédito de injeção (Tabela 7 e fronteira)
        'contrato_sweep_kw': list(range(400, 801, 50)),   # eixo contrato de demanda (Tabela 9)
        'verificar_fronteira_milp': True,
        'pasta_saida': '.',
        'albedo': 0.20, 'percentil_critico': 10.0, 'serie_fv_origem': 'NASA_POWER',
        'seed_base': 42, 'retomar': False,
        'seed_solver': 42, 'threads_solver': 1, 'ano_financeiro': 2027,
        'feriados_sem_ponta': [], 'capacidades_bess_kwh': [500., 750., 1000.],
        'fom_bess_fracao_ano': 0.025, 'limite_total_caso_s': 900,
        'alvo_exato': True,
        'tol_potencia_kw': 1e-3, 'tol_energia_kwh': 1e-3, 'tol_soc': 1e-6,
        'tol_fatura_rs': 0.05,
        'referencia_temporal': 'LST NASA; aproximação hora solar para geometria; ponta 18–21 na mesma grade',

    }
    return p


def _q_perm(p, m_h):
    q = dict(p)
    q.update({'perm_media_h': m_h, 'perm_min_h': max(1.0, m_h - 4.0), 'perm_max_h': min(20.0, m_h + 4.0),
              'perm_desvio_h': min(p['perm_desvio_h'], 2.0), 'limitar_partida_manha': False})
    return q


def _jsonavel(x):
    if isinstance(x, dict):
        return {str(k): _jsonavel(v) for k, v in x.items()}
    if isinstance(x, (tuple, list, np.ndarray)):
        return [_jsonavel(v) for v in x]
    if isinstance(x, (np.integer, np.bool_)):
        return x.item()
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    return x


def _hash_obj(obj):
    return hashlib.sha256(json.dumps(_jsonavel(obj), sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def _hash_codigo():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def versoes_ambiente():
    out = {'python': platform.python_version(), 'plataforma': platform.platform()}
    for pacote in ('pulp', 'highspy', 'numpy', 'pandas', 'matplotlib'):
        try:
            out[pacote] = version(pacote)
        except PackageNotFoundError:
            out[pacote] = 'não instalado'
    return out


def _status_modelo(m, ret):
    # https://coin-or.github.io/pulp/technical/constants.html
    status = pulp.LpStatus.get(m.status, 'Undefined')
    sol = getattr(m, 'sol_status', None)
    tem = sol in (pulp.LpSolutionOptimal, pulp.LpSolutionIntegerFeasible)
    return status, tem


def resolver_lcos(p):
    # Não converter os antigos 6.000 ciclos sem fonte em indicador científico.
    return float('nan')


def _tabela_fv(p):
    df, meta = carregar_solar_completo(p)
    grupos = ((12, 1, 2), (6, 7, 8))
    tab = np.array([df[df.MO.isin(meses)].groupby('HR').POA.mean().reindex(range(24)).to_numpy()
                    for meses in grupos])
    if not np.isfinite(tab).all():
        raise ValueError('A série solar não cobre todas as horas das duas estações.')
    return tab, f"NASA POWER 2019–2023 (período configurado: {p['anos_serie_fv']}), LST; fonte {meta['hash_csv']}"


def carregar_solar_completo(p):
    # NASA POWER: https://power.larc.nasa.gov/docs/services/api/temporal/hourly/
    # LST NASA é uma faixa de longitude de 15 graus, não hora solar aparente exata.
    # Mantém-se a aproximação temporal declarada no E4; a série é horária, não medida a 15 min.
    parametros = {k: p[k] for k in ('lat_graus', 'lon_graus', 'anos_serie_fv',
                                   'inclinacao_fv_graus', 'albedo')}
    parametros.update(fonte='NASA POWER hourly', time_standard='LST',
                      variaveis=['GHI', 'DNI', 'DIF'], transposicao='isotropica_V21_E4')
    chave = _hash_obj(parametros)
    if chave in _CACHE_SOLAR_COMPLETO:
        return _CACHE_SOLAR_COMPLETO[chave]
    pasta = Path(p['arquivo_serie_fv']).expanduser().resolve()
    if pasta.suffix.lower() == '.csv':
        raise ValueError('V22 exige série horária completa. Configure arquivo_serie_fv como pasta de dados, não o CSV sazonal antigo.')
    pasta.mkdir(parents=True, exist_ok=True)
    arquivo, manifesto = pasta/(chave+'.csv'), pasta/(chave+'.json')
    if arquivo.exists() and manifesto.exists():
        meta = json.loads(manifesto.read_text(encoding='utf8'))
        if meta['parametros'] != parametros or meta['hash_csv'] != hashlib.sha256(arquivo.read_bytes()).hexdigest():
            raise ValueError('Cache solar divergente ou modificado; revise os arquivos e a origem.')
        df = pd.read_csv(arquivo)
    else:
        if not p.get('baixar_serie_fv', True):
            raise FileNotFoundError('Série horária NASA ausente. Habilite download ou disponibilize o cache completo com manifesto.')
        print('  Obtendo GHI, DNI e difusa da NASA POWER. Não haverá substituição por curva sintética.')
        df = baixar_nasa_power_horario(p['lat_graus'], p['lon_graus'], p['anos_serie_fv'], True)
        tmp_solar = arquivo.with_suffix('.csv.tmp')
        df.to_csv(tmp_solar, index=False)
        substituir_com_retry(tmp_solar, arquivo)
        meta = dict(parametros=parametros, coletado_utc=datetime.now(timezone.utc).isoformat(),
                    hash_csv=hashlib.sha256(arquivo.read_bytes()).hexdigest(),
                    unidades='GHI/DNI/DIF: W/m² (média horária)',
                    url='https://power.larc.nasa.gov/api/temporal/hourly/point')
        tmp_meta = manifesto.with_suffix('.json.tmp')
        tmp_meta.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf8')
        substituir_com_retry(tmp_meta, manifesto)
    df['data'] = pd.to_datetime(dict(year=df.YEAR, month=df.MO, day=df.DY), errors='coerce')
    validas, descartadas = [], []
    for data, g in df.groupby('data', dropna=False):
        a = g[['GHI','DNI','DIF']].to_numpy(float)
        ok = (not pd.isna(data) and len(g)==24 and set(g.HR)==set(range(24))
              and np.isfinite(a).all() and (a >= 0).all())
        if ok:
            validas.append(g)
        else:
            descartadas.append({'data': str(data), 'linhas': len(g), 'motivo': 'dia incompleto, duplicado ou irradiância inválida'})
    if not validas:
        raise ValueError('Nenhum dia solar completo e válido.')
    df = pd.concat(validas).sort_values(['data','HR']).reset_index(drop=True)
    presentes=set(df['data'].dt.strftime('%Y-%m-%d'))
    ja_descartados={x['data'][:10] for x in descartadas}
    for ano in p['anos_serie_fv']:
        for data in pd.date_range(f'{ano}-01-01',f'{ano}-12-31'):
            dia=data.strftime('%Y-%m-%d')
            if dia not in presentes and dia not in ja_descartados:
                descartadas.append({'data':dia,'linhas':0,'motivo':'dia inteiro ausente'})
    if set(df.YEAR.astype(int)) != set(p['anos_serie_fv']):
        raise ValueError('Falta pelo menos um ano solicitado na série solar válida.')
    df['POA'] = irradiancia_plano_inclinado(df, p['lat_graus'], p['inclinacao_fv_graus'], p['albedo'])
    diario = df.groupby('data', as_index=False).POA.sum().rename(columns={'POA':'Hi_kwh_m2_dia'})
    diario.Hi_kwh_m2_dia /= 1000.0
    meta.update(dias_validos=len(diario), dias_descartados=len(descartadas),
                Hi_anual_media=float(diario.Hi_kwh_m2_dia.mean()),
                atlas_comparacao='Sudeste, plano inclinado, 5,26 kWh/m².dia; Atlas 2017 p.67 Fig.65; não calibra o ponto local')
    salvar_csv(diario, 'solar_energia_diaria.csv')
    if descartadas:
        salvar_csv(pd.DataFrame(descartadas, columns=['data','linhas','motivo']), 'solar_dias_descartados.csv')
    meta['arquivo_cache'] = str(arquivo.resolve())
    gravar_json('solar_proveniencia.json', meta)
    _CACHE_SOLAR_COMPLETO[chave] = df, meta
    return df, meta


def auditar_solucao(r, frota, p, rotulo, cap_bess):
    """Auditoria independente, por visita e por intervalo. Nunca altera o status do solver."""
    av = []
    tk, te, ts = p['tol_potencia_kw'], p['tol_energia_kwh'], p['tol_soc']
    P = np.asarray(r['P_veic'])
    campos = [P, r['imp'], r['exp'], r['pv_gen'], r['bess_in'], r['bess_out'], r['bess_soc']]
    if not all(np.isfinite(x).all() for x in campos):
        r['auditoria_ok'] = False
        r['auditoria_avisos'] = ['ausência de solução numérica completa']
        return r['auditoria_avisos']
    janela = taper = excesso_soc = 0.0
    finais = []
    for v in range(frota['f']):
        ids = (int(frota['a'][v]) + np.arange(int(frota['L'][v]))) % T_STEPS
        fora = np.ones(T_STEPS, bool); fora[ids] = False
        if fora.any():
            janela = max(janela, float(np.max(np.abs(P[v, fora]))))
        soc = float(frota['soc_i'][v])
        for t in ids:
            taper = max(taper, float(P[v,t] - pmax_carregador(soc, p)))
            soc += p['eta_carga'] * DT * P[v,t] / frota['cap'][v]
            excesso_soc = max(excesso_soc, soc-1, -soc)
        finais.append(soc)
    finais = np.array(finais)
    deficit = float(np.maximum(0, (p['soc_f']-finais)*frota['cap']).sum())
    excesso_alvo = float(np.maximum(0, (finais-p['soc_f'])*frota['cap']).sum())
    rec = P.sum(axis=0)
    bal = r['imp']+r['pv_gen']+r['bess_out']-rec-r['exp']-r['bess_in']
    e_dinamica = (finais-frota['soc_i'])*frota['cap']
    e_rec_bat = P.sum(axis=1)*DT*p['eta_carga']
    erro_veic = float(np.max(np.abs(e_rec_bat-e_dinamica)))
    erro_alvo_ca = float(rec.sum()*DT - frota['e_nec'].sum()/p['eta_carga'])
    dtb = p['eta_bess']*r['bess_in']*DT - r['bess_out']*DT/p['eta_bess']
    energia_b = np.asarray(r['bess_soc'])*cap_bess/100.0
    erro_b_din = float(np.max(np.abs(energia_b-np.roll(energia_b,1)-dtb))) if cap_bess>0 else 0.0
    ciclo = float(dtb.sum())
    n = r.get('N_aud', p['n_carregadores'])
    simult = int(np.max((P>tk).sum(axis=0)))
    relaxado = r.get('ocupacao_relaxada',False)
    verificacoes = {
        'potência negativa': min(float(np.min(x)) for x in (P,r['imp'],r['exp'],r['bess_in'],r['bess_out'])) < -tk,
        'potência fora da janela': janela>tk,
        'limite/taper do carregador': taper>tk or float(P.max())>p['pot_carregador']+tk,
        'SoC veicular fora de [0,1]': excesso_soc>ts,
        'SoC alvo não atendido': deficit>te,
        'energia acima do alvo exato': excesso_alvo>te,
        'simultaneidade inteira': not relaxado and simult>n,
        'capacidade agregada dos carregadores': float(rec.max())>n*p['pot_carregador']+tk,
        'balanço de potência': float(np.max(np.abs(bal)))>tk,
        'balanço de energia no ciclo': abs(float(bal.sum()*DT))>te,
        'dinâmica veicular': erro_veic>te,
        'dinâmica do BESS': erro_b_din>te,
        'ciclo do BESS': abs(ciclo)>te,
        'carga e descarga simultâneas': float(np.minimum(r['bess_in'],r['bess_out']).max())>tk,
        'importação e exportação simultâneas': float(np.minimum(r['imp'],r['exp']).max())>tk,
        'FIFO não periódico': r.get('residuo_periodico_kw',0)>tk,
    }
    if cap_bess>0:
        verificacoes['limites BESS'] = (np.min(r['bess_soc'])/100 < p['soc_min_bess']-ts or
            np.max(r['bess_soc'])/100 > p['soc_max_bess']+ts or
            max(np.max(r['bess_in']),np.max(r['bess_out'])) > p['pot_bess']+tk)
    else:
        verificacoes['BESS ausente com fluxo'] = max(np.max(r['bess_in']),np.max(r['bess_out']))>tk
    av = [k for k,v in verificacoes.items() if v]
    residuos = dict(janela_kw=janela,taper_kw=taper,balanco_kw=float(np.max(np.abs(bal))),
                    balanco_kwh=float(bal.sum()*DT),dinamica_veiculo_kwh=erro_veic,
                    energia_alvo_ca_residuo_kwh=erro_alvo_ca,deficit_kwh=deficit,
                    excesso_alvo_kwh=excesso_alvo,bess_dinamica_kwh=erro_b_din,bess_ciclo_kwh=ciclo,
                    simultaneidade_max=simult,ocupacao_relaxada=relaxado)
    r.update(auditoria_ok=not av,auditoria_avisos=av,auditoria=residuos,
             energia_nao_atendida=deficit,n_nao_atendidos=int(np.sum((p['soc_f']-finais)*frota['cap']>te)),
             soc_final=finais,cronograma_fisico=(not av and simult<=n))
    REGISTRO_AUDITORIAS.append(dict(rotulo=rotulo,**residuos,aprovado=not av,avisos=' | '.join(av)))
    if av:
        print(f"  [AUDITORIA] {rotulo}: " + '; '.join(av))
    return av


def executar_verificacao(p):
    print('\nVERIFICAÇÃO A–F — parâmetros analíticos fixos do E4; não dependem das escolhas interativas.')
    q=parametros_padrao()
    q.update(pot_carregador=50.0,n_carregadores=1,cap_bat=100.0,soc_f=.9,usar_ccv=False,
             eta_carga=.95,contrato_modo='otimo',modo_simultaneidade='milp',solver=p['solver'],
             tempo_limite_s=p['tempo_limite_s'],tempo_max_s=p['tempo_max_s'],gap_rel=p['gap_rel'])
    def frota(n,L,a=88):
        return dict(f=n,a=np.full(n,a),L=np.full(n,L),cap=np.full(n,100.),soc_i=np.full(n,.2),e_nec=np.full(n,70.))
    casos=[('A',frota(1,16),q,0,None),('B',frota(1,4),q,0,None),('C',frota(2,16),q,0,None),
           ('D',frota(1,4),dict(q,pot_carregador=100.,pot_bess=50.),100.,None),
           ('E',frota(1,16),q,0,10.),('F',frota(1,12,72),dict(q,pot_carregador=100.),0,None)]
    linhas=[]
    def medir(caso,grandeza,esperado,obtido,tol,certificado):
        ok=bool(certificado and np.isfinite(obtido) and abs(obtido-esperado)<=tol)
        linhas.append(dict(caso=caso,grandeza=grandeza,esperado=esperado,obtido=obtido,
                           erro=obtido-esperado,tolerancia=tol,aprovado=ok))
        if not ok:
            print(f'  FALHA {caso}: {grandeza}; esperado {esperado:.6f}, obtido {obtido:.6f}.')
    for cod,f,qp,cb,D in casos:
        r=resolver_ems(f,qp,cb,0,np.zeros(24),'verde','milp',1,D,'Verif-'+cod)
        tem=r.get('tem_solucao',False)
        if cod=='B':
            medir(cod,'déficit (kWh)',22.5,r.get('energia_nao_atendida',np.nan),.001,
                  tem and r.get('infeasivel_original') and r.get('convergido'))
        else:
            res=_aval_otim(r,cb,qp,'verde',D,f)
            cert=res.get('valido',False) and res.get('convergido',False)
            if cod in ('A','C','D'):
                esperado={'A':70/.95/4,'C':140/.95/4,'D':70/.95-50}[cod]
                medir(cod,'pico (kW)',esperado,res['pico'],.05,cert)
            if cod=='D':
                medir(cod,'descarga BESS (kWh)',50,float(np.sum(res['bess_out'])*DT),.05,cert)
                medir(cod,'reposição BESS (kWh)',50/.95**2,float(np.sum(res['bess_in'])*DT),.05,cert)
            if cod=='E':
                medir(cod,'demanda (R$/mês)',12.02*70/.95/4,res['custo_demanda'],.05,cert)
                medir(cod,'ultrapassagem (R$/mês)',2*12.02*(70/.95/4-10),res['custo_ultrapassagem'],.05,cert)
                medir(cod,'demanda + ultrapassagem (R$/mês)',12.02*70/.95/4+2*12.02*(70/.95/4-10),res['custo_demanda']+res['custo_ultrapassagem'],.05,cert)
            if cod=='F':
                medir(cod,'energia (R$/mês)',(70/.95)*(22*1.81466+8*.45128),res['custo_energia'],.05,cert)
        salvar_csv(pd.DataFrame(linhas),'verificacao_analitica.csv')
    # Casos suplementares: descontinuidade de faturamento exatamente na banda e acima.
    for rot,pico in (('E-limite',10*1.05),('E-acima',10*1.05+.001)):
        dem,ult=custo_demanda_verde_fixo(pico,10,q)
        medir(rot,'ultrapassagem (R$/mês)',0 if rot=='E-limite' else 2*12.02*(pico-10),ult,1e-8,True)
    salvar_csv(pd.DataFrame(linhas),'verificacao_analitica.csv')
    aprovado=all(l['aprovado'] for l in linhas)
    print(f'Verificação: {sum(l["aprovado"] for l in linhas)}/{len(linhas)} verificações aprovadas. Detalhes em tabelas/verificacao.csv.')
    gravar_json('verificacao_status.json',dict(aprovado=aprovado,codigo=_hash_codigo(),
                                              verificacoes=len(linhas),rodada_utc=datetime.now(timezone.utc).isoformat()))
    return aprovado

# -----------------------------------------------------------------------------
# Saída enxuta: sobrescritas atômicas de tabelas consolidadas, sem arquivos por caso.
# -----------------------------------------------------------------------------
def caminho(nome):
    dest=Path(_PASTA['saida'])/nome
    dest.parent.mkdir(parents=True,exist_ok=True)
    return dest


def substituir_com_retry(tmp,dest):
    for tentativa in range(12):
        try:
            os.replace(tmp,dest)
            return
        except PermissionError:
            if tentativa==11:
                raise PermissionError(f'Feche o arquivo {dest.name} no Excel/visualizador. Temporário preservado: {tmp}')
            time.sleep(.25*(tentativa+1))


def gravar_json(nome,obj):
    if nome in ('solar_proveniencia.json','verificacao_status.json'):
        REGISTRO[nome.removesuffix('.json')]=copy.deepcopy(obj)
        return
    dest=caminho(nome);tmp=dest.with_suffix(dest.suffix+'.tmp')
    tmp.write_text(json.dumps(_jsonavel(obj),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
    substituir_com_retry(tmp,dest)


def salvar_csv(df,nome):
    if nome=='solar_energia_diaria.csv':
        return  # série integral já no cache de entrada, referenciado por hash
    if nome=='solar_dias_descartados.csv':
        REGISTRO['dias_solares_descartados']=df.to_dict('records')
        return
    if nome=='verificacao_analitica.csv':
        VERIFICACOES.clear();VERIFICACOES.extend(df.to_dict('records'))
        nome='tabelas/verificacao.csv'
    dest=caminho(nome);tmp=dest.with_suffix(dest.suffix+'.tmp')
    df.to_csv(tmp,index=False,sep=';',decimal=',',na_rep='N/D',encoding='utf-8-sig')
    substituir_com_retry(tmp,dest)


def numero(x):
    if not _finito(x):return 'N/D'
    if abs(x)<.005:x=0.
    return f'{x:,.2f}'.replace(',','@').replace('.',',').replace('@','.')


def publicar():
    for rows,nome in ((LINHAS,'cenarios'),(COMPARACOES,'comparacoes'),
                      (FINANCEIRO,'financeiro_anual'),(FRONTEIRA,'fronteira'),
                      (DIMENSIONAMENTO,'capacidade_bess'),(VERIFICACOES,'verificacao')):
        if rows:salvar_csv(pd.DataFrame(rows),'tabelas/'+nome+'.csv')
    if AUDITORIAS:salvar_csv(pd.DataFrame(AUDITORIAS),'auditoria/auditoria.csv')
    if REGISTRO_SOLVES:salvar_csv(pd.DataFrame(REGISTRO_SOLVES),'auditoria/solver.csv')
    if PERFIS:salvar_csv(pd.DataFrame(PERFIS),'tabelas/perfis_operacao.csv')
    gravar_json('auditoria/registro_execucao.json',REGISTRO)


class HiGHSComInicio(pulp.HiGHS):
    """PuLP 3.3.2 callSolver hook; highspy 1.15.1 setSolution sparse overload.
    Código das interfaces consultado estaticamente; não testado pelo editor.
    kOk confirma recebimento pela API, não viabilidade nem uso do incumbente.
    """
    def __init__(self,*args,enviar_inicio=False,**kwargs):
        super().__init__(*args,**kwargs)
        self.enviar_inicio=enviar_inicio
        self.inicio_enviado=False
        self.inicio_status='não solicitado'

    def callSolver(self,lp):
        if self.enviar_inicio:
            try:
                pares=[(v.index,float(v.varValue)) for v in lp.variables()
                       if v.varValue is not None and np.isfinite(v.varValue)]
                if pares:
                    ids=np.array([a for a,b in pares],dtype=np.int32)
                    vals=np.array([b for a,b in pares],dtype=np.float64)
                    self.inicio_enviado=True
                    status=lp.solverModel.setSolution(len(pares),ids,vals)
                    self.inicio_status=str(status)+'; recebimento pela API, viabilidade não certificada'
            except (AttributeError,TypeError,ValueError,RuntimeError) as exc:
                self.inicio_status='não aplicado: '+str(exc)
        super().callSolver(lp)


def criar_solver(p,tempo_limite=None,warm=False):
    if p['solver']!='HiGHS':
        raise ValueError('Esta versão usa HiGHS/highspy; não troca o solver silenciosamente.')
    solver=HiGHSComInicio(msg=bool(p['solver_msg']),timeLimit=tempo_limite or p['tempo_limite_s'],
                         gapRel=p['gap_rel'],threads=p['threads_solver'],
                         random_seed=p['seed_solver'],enviar_inicio=warm)
    if not solver.available():raise RuntimeError('HiGHS indisponível; instale highspy no mesmo Python.')
    return solver,'HiGHS'


def _auditar_raw(r,frota,p,cb):
    if not r.get('tem_solucao'):return False
    z=dict(r,imp=r['grid_import'],exp=r['grid_export'])
    auditar_solucao(z,frota,p,'candidato',cb)
    r['auditoria_ok']=z['auditoria_ok']
    r['auditoria_avisos']=z['auditoria_avisos']
    r['auditoria']=z['auditoria']
    return z['auditoria_ok']


def resolver_ems(frota,p,cap_bess,pot_fv,curva,modalidade='verde',modo=None,
                 n_carregadores=None,dem_contratada=None,rotulo='',warm=None):
    _ID_CHAMADA[0]+=1;cid=_ID_CHAMADA[0];inicio=time.perf_counter()
    candidatos=[];outros=[];bounds=[];ultima={};tent=0
    limite=min(p['tempo_limite_s'],p['tempo_max_s'])
    while True:
        restante=p['limite_total_caso_s']-(time.perf_counter()-inicio)
        if restante<1 and tent:break
        tent+=1;tl=max(1.,min(limite,restante))
        r=otimizar_ems(frota,p,cap_bess,pot_fv,curva,modalidade,modo,n_carregadores,
                      dem_contratada,rotulo,warm,tempo_limite=tl,tentativa=tent,id_chamada=cid)
        ultima=r
        REGISTRO_SOLVES[-1]['bloco']=CONTEXTO_ATUAL or rotulo
        if _finito(r.get('bound')):bounds.append(r['bound'])
        if r.get('tem_solucao'):
            outros.append(r)
            if _auditar_raw(r,frota,p,cap_bess):
                candidatos.append(r)
                warm=dict(r,bess_in=r['bess_in'],bess_out=r['bess_out'])
        if r.get('infeasivel'):
            if candidatos:
                r=dict(min(candidatos,key=lambda x:x['obj']),conflito_solver=True,convergido=False)
                r['status']='Certificados incompatíveis';return r
            diag=otimizar_ems(frota,p,cap_bess,pot_fv,curva,modalidade,modo,n_carregadores,
                             dem_contratada,rotulo+' [diagnóstico]',permitir_deficit=True,
                             tempo_limite=p['tempo_max_s'],tentativa=tent+1,id_chamada=cid)
            REGISTRO_SOLVES[-1]['bloco']=CONTEXTO_ATUAL or rotulo
            diag.update(infeasivel_original=True,status_estrito='Infeasible',
                        tempo_total_s=time.perf_counter()-inicio,tentativas_total=tent+1)
            return diag
        if r.get('convergido') and r.get('auditoria_ok'):break
        if not p['escalonar_tempo'] or limite>=p['tempo_max_s']:break
        limite=min(p['tempo_max_s'],limite*3)
    r=copy.deepcopy(min(candidatos,key=lambda x:x['obj']) if candidatos else
                    min(outros,key=lambda x:x['obj']) if outros else ultima)
    if candidatos and bounds:
        lb=max(bounds);ub=r['obj']
        if lb>ub+p['tol_fatura_rs']:
            r.update(conflito_solver=True,convergido=False,status='Limites incompatíveis')
        else:
            r['bound']=min(lb,ub)
            r['gap']=max(0.,ub-lb)/max(abs(ub),1e-10)
            r['convergido']=r['gap']<=p['gap_rel']+1e-9
    r.update(tempo_total_s=time.perf_counter()-inicio,tentativas_total=tent,infeasivel_original=False)
    return r


def _aval_otim(r,cb,p,modalidade,D,frota):
    if not r.get('tem_solucao'):
        z=np.full(96,np.nan)
        out=dict(r,imp=z,exp=z,pv_gen=z,recarga=z,bess_in=z,bess_out=z,bess_soc=z,
                 P_veic=np.full((frota['f'],96),np.nan),valido=False,auditoria_ok=False,
                 opex=np.nan,pico=np.nan,energia_nao_atendida=np.nan,
                 auditoria_avisos=['sem solução utilizável'])
        return out
    out=avaliar(r['grid_import'],r['grid_export'],r['pv_gen'],r['bess_out'],cb,p,
                modalidade,D,extras={k:v for k,v in r.items() if k not in ('grid_import','grid_export')})
    auditar_solucao(out,frota,p,CONTEXTO_ATUAL,cb)
    erro=abs(r['obj']-out['opex']-out['regularizador_rs_mes']) if not r.get('diagnostico') else np.nan
    out['residuo_objetivo_rs']=erro
    if _finito(erro) and erro>p['tol_fatura_rs']:
        out['auditoria_ok']=False;out['auditoria_avisos'].append('objetivo não reconciliado com fatura e regularizador')
    out['valido']=bool(out['auditoria_ok'] and not r.get('diagnostico') and
                       not r.get('infeasivel_original') and not r.get('conflito_solver'))
    out['opex_bruto']=out['opex']
    if not out['valido']:out['opex']=np.nan
    out['atendimento_pct']=100*(1-out['energia_nao_atendida']/max(frota['e_nec'].sum(),1e-12))
    # Limites sobre o MENOR custo de fatura do conjunto viável, não sobre a fatura
    # do ótimo de J=B+R. L_J-R_max <= min B <= B(incumbente). Sem diagnóstico.
    rmax=_c_deg(p)*p['pot_bess']*24*p['dias_mes'] if cb else 0.
    out['fatura_min_LB']=max(0.,r['bound']-rmax) if out['valido'] and _finito(r.get('bound')) else np.nan
    out['fatura_min_UB']=out['opex'] if out['valido'] else np.nan
    out['limite_fatura_tipo']='mínimo de B; conservador, não ótimo de B+R' if cb else 'mesma função objetivo e fatura'
    return out


def estado(r):
    if r.get('valido'):
        return 'HEURISTICA_VIAVEL' if r.get('status')=='Heuristica' else ('CONVERGIDO' if r.get('convergido') else 'VIAVEL_COM_GAP')
    if r.get('infeasivel_original'):return 'INVIAVEL_COMPROVADO'
    if r.get('status')=='Heuristica':return 'HEURISTICA_COM_DEFICIT'
    if not r.get('tem_solucao'):return 'SEM_SOLUCAO'
    return 'AUDITORIA_REPROVADA'


def linha(r,p,contexto):
    get=lambda k:r.get(k,np.nan)
    ok=r.get('valido',False)
    return dict(bloco=contexto,cenario=r['cenario'],modalidade=r['modalidade'],seed=r['seed'],
        frota=r['f'],carregadores=r['N'],permanencia_nominal_h=p['perm_media_h'],
        permanencia_efetiva_h=get('permanencia_efetiva_h'),credito=p['multiplicador_credito'],
        capacidade_BESS_kWh=p['cap_bess'] if r['cenario']=='C3' else 0,
        contrato_kW=r.get('contrato_kw'),estado=estado(r),status_solver=get('status'),
        gap_pct=100*get('gap'),tempo_total_s=get('tempo_total_s'),tentativas=get('tentativas_total'),
        objetivo_LB=get('bound'),objetivo_UB=get('obj'),
        fatura_min_LB=get('fatura_min_LB'),fatura_min_UB=get('fatura_min_UB'),
        limite_fatura_tipo=r.get('limite_fatura_tipo','N/A'),
        pico_kW=get('pico') if ok else np.nan,fator_carga=get('fator_carga') if ok else np.nan,
        fatura_Rmes=get('opex'),demanda_Rmes=get('custo_demanda') if ok else np.nan,
        ultrapassagem_Rmes=get('custo_ultrapassagem') if ok else np.nan,
        energia_Rmes=get('custo_energia') if ok else np.nan,regularizador_Rmes=get('regularizador_rs_mes') if ok else np.nan,
        autoconsumo_FV_pct=get('autoconsumo') if ok else np.nan,injecao_kWh_dia=get('injecao_kwh_dia') if ok else np.nan,
        EFC_mes=get('efc') if ok else np.nan,energia_ponta_kWh_dia=get('ponta_kwh') if ok else np.nan,
        deficit_kWh_dia=get('energia_nao_atendida'),atendimento_pct=get('atendimento_pct'),
        espera_C0_min=60*get('espera_media_h') if r['cenario'] in ('C0','C0b') else np.nan,
        ocupacao_relaxada=bool(r.get('ocupacao_relaxada')),aceito_fisicamente=ok,
        avisos='; '.join(r.get('auditoria_avisos',[])),chave=r['chave'])


def grupo(contexto,p,curva,f=None,N=None,seed=None,modo='milp',D=None,
          modalidades=('verde',),codigos=('C0','C0b','C1','C2','C3'),perfil=False):
    global CONTEXTO_ATUAL
    CONTEXTO_ATUAL=contexto
    seed=p['seed_base'] if seed is None else seed
    print(f'\nCalculando {contexto} ({modo.upper()})...',flush=True)
    f=int(f or p['frota_base']);N=int(N or p['n_carregadores']);frota=gerar_frota(f,p,seed)
    out={}
    for mod in modalidades:
        for cod in codigos:
            ident=dict(parametros={k:v for k,v in p.items() if k not in ('pasta_saida','retomar')},
                       curva=curva,f=f,N=N,seed=seed,modo=modo,modalidade=mod,cenario=cod,D=D)
            chave=_hash_obj(ident)
            if chave in MEMORIA:r=copy.deepcopy(MEMORIA[chave])
            else:
                inicio=time.perf_counter();cb=p['cap_bess'] if cod=='C3' else 0.
                if cod in ('C0','C0b'):
                    s=simular_nao_coordenada(frota,p,N,cod=='C0b');z=np.zeros(96)
                    r=avaliar(s['carga'],z,z,z,0,p,mod,D,
                              extras=dict(s,recarga=s['carga'],bess_in=z,bess_out=z,bess_soc=z,pv_gen=z,
                                          N_aud=N,ocupacao_relaxada=False,status='Heuristica',
                                          convergido=False,tem_solucao=True))
                    auditar_solucao(r,frota,p,contexto,0)
                    r['valido']=r['auditoria_ok'];r['opex_bruto']=r['opex']
                    if not r['valido']:r['opex']=np.nan
                    r['atendimento_pct']=100*(1-r['energia_nao_atendida']/frota['e_nec'].sum())
                else:
                    antes={'C1':'C0','C2':'C1','C3':'C2'}[cod]
                    warm=out.get((mod,antes))
                    if warm and not warm.get('valido'):warm=None
                    raw=resolver_ems(frota,p,cb,p['pot_fv'] if cod in ('C2','C3') else 0.,curva,
                                     mod,modo,N,D,f'{contexto}/{mod}/{cod}',warm)
                    r=_aval_otim(raw,cb,p,mod,D,frota)
                r.setdefault('tempo_total_s',time.perf_counter()-inicio)
                r.update(cenario=cod,modalidade=mod,seed=seed,f=f,N=N,chave=chave,contrato_kw=D,
                         permanencia_efetiva_h=float(frota['perm_h'].mean()))
                MEMORIA[chave]=copy.deepcopy(r)
            out[mod,cod]=r
            row=linha(r,p,contexto)
            ref=out.get((mod,'C0'))
            row['reducao_pico_vs_C0_pct']=(100*(ref['pico']-r['pico'])/ref['pico']
                if ref and ref.get('valido') and r.get('valido') and ref['pico']>0 else np.nan)
            row['D_sem_ultrapassagem_para_pico_kW']=(r['pico']/(1+p['tolerancia_ultrapassagem'])
                if r.get('valido') and mod=='verde' else np.nan)
            row['energia_BESS_descarga_kWh_dia']=float(np.sum(r['bess_out'])*DT) if r.get('valido') else np.nan
            LINHAS.append(row)
            a=r.get('auditoria',{})
            AUDITORIAS.append(dict(bloco=contexto,modalidade=mod,cenario=cod,chave=chave,
                estado=estado(r),**a,
                energia_requerida_bateria_kWh=float(frota['e_nec'].sum()),
                energia_entregue_bateria_kWh=float(np.sum(r['P_veic'])*DT*p['eta_carga']),
                E_FV_kWh=float(np.sum(r['pv_gen'])*DT),E_imp_kWh=float(np.sum(r['imp'])*DT),
                E_exp_kWh=float(np.sum(r['exp'])*DT),E_ch_kWh=float(np.sum(r['bess_in'])*DT),
                E_dis_kWh=float(np.sum(r['bess_out'])*DT),N_carregadores=N))
            if perfil and r.get('valido'):
                for t in range(96):
                    PERFIS.append(dict(bloco=contexto,modalidade=mod,cenario=cod,hora=t*DT,
                        rede_kW=r['imp'][t],exportacao_kW=r['exp'][t],FV_kW=r['pv_gen'][t],
                        recarga_kW=r['recarga'][t],BESS_carga_kW=r['bess_in'][t],
                        BESS_descarga_kW=r['bess_out'][t],BESS_SoC_pct=r['bess_soc'][t]))
    for mod in modalidades:
        for ca,cb,tipo in (('C0','C1','gestao'),('C0','C0b','regra_simples'),('C0b','C1','otimizacao'),
                           ('C1','C2','FV'),('C2','C3','BESS')):
            if (mod,ca) not in out or (mod,cb) not in out:continue
            a,b=out[mod,ca],out[mod,cb];d=diferenca(a,b)
            COMPARACOES.append(dict(bloco=contexto,modalidade=mod,comparacao=tipo,
                economia_Rmes=d,reducao_pico_kW=a.get('pico',np.nan)-b.get('pico',np.nan) if _finito(d) else np.nan,
                estados=estado(a)+' / '+estado(b),
                observacao='estimativa de incumbentes; consultar limites' if _finito(d) else 'não comparável',
                diferenca_negativa=bool(_finito(d) and d < -p['tol_fatura_rs'])))
    print(f'  {contexto}: '+', '.join(f'{m}/{c} {estado(r)}' for (m,c),r in out.items()))
    REGISTRO['ultimo_bloco_concluido']=contexto
    publicar()
    return out


def diferenca(a,b):
    return a['opex']-b['opex'] if a.get('valido') and b.get('valido') else np.nan


def preparar_entradas(p):
    df,meta=carregar_solar_completo(p)
    tab,_=_tabela_fv(p);curvas={'verao':tab[0]/1000,'inverno':tab[1]/1000}
    inverno=df[df.MO.isin([6,7,8])]
    diaria=inverno.groupby('data').POA.sum().sort_index()/1000
    p5,p10,p15=np.percentile(diaria,[5,10,15],method='linear')
    dias=diaria.index[(diaria>=p5)&(diaria<=p15)]
    if not len(dias):raise ValueError('Sem dias completos P5–P15.')
    observado=(diaria-p10).abs().idxmin()  # empate: primeira data na ordenação crescente
    curvas['p5_p15']=inverno[inverno.data.isin(dias)].groupby('HR').POA.mean().reindex(range(24)).to_numpy()/1000
    curvas['p10_observado']=inverno[inverno.data==observado].sort_values('HR').POA.to_numpy()/1000
    if any(len(c)!=24 or not np.isfinite(c).all() for c in curvas.values()):raise ValueError('Perfil solar incompleto.')
    REGISTRO['perfis_solares']=dict(fonte='NASA_POWER',ano_inicio=min(p['anos_serie_fv']),ano_fim=max(p['anos_serie_fv']),
        P5=p5,P10=p10,P15=p15,criterio='energia diária POA de dias completos de junho–agosto',
        dias_media=[str(x.date()) for x in dias],dia_observado=str(observado.date()),
        energia_observado_kWh_m2=float(diaria.loc[observado]),distancia_P10=float(abs(diaria.loc[observado]-p10)),
        empate='primeira data cronológica',replicacao='cada hora repetida em quatro intervalos de 15 min',
        limite='NASA horária suaviza rampas; estes casos não medem a máxima demanda mensal real')
    REGISTRO['atlas']=dict(Hi_NASA=float(df.groupby('data').POA.sum().mean()/1000),
        H_NASA=float(df.groupby('data').GHI.sum().mean()/1000),Hi_regional_INPE=5.26,H_regional_INPE=5.06,
        origem='Atlas 2017, p.67 Fig.65; plano inclinado vs horizontal; não calibração local')
    salvar_csv(pd.DataFrame([dict(perfil=k,hora=h,POA_Wm2=v*1000) for k,c in curvas.items()
                            for h,v in enumerate(c)]),'entrada/perfis_solares.csv')
    f=gerar_frota(p['frota_base'],p,p['seed_base']);g=gerar_frota(100,p,p['seed_base'])
    REGISTRO['frota_prefixo_54_100']=bool(all(np.array_equal(f[k],g[k][:p['frota_base']])
                                                   for k in ('a','L','soc_i','cap')))
    salvar_csv(pd.DataFrame(dict(veiculo=range(1,f['f']+1),chegada_h=f['chegada_h'],
        chegada_grade_h=f['a']*DT,partida_h=f['partida_h'],permanencia_h=f['perm_h'],
        SoC_inicial=f['soc_i'],energia_bateria_kWh=f['e_nec'])),'entrada/frota_base.csv')
    return curvas


def calendario_anual(p):
    # Hipótese explícita: out–mar usa perfil DJF; abr–set usa perfil JJA.
    # Calendário de dias reais; feriados informados pelo usuário são retirados da ponta.
    feriados={date.fromisoformat(x) for x in p['feriados_sem_ponta']}
    out={}
    for est,meses in [('verao',[1,2,3,10,11,12]),('inverno',[4,5,6,7,8,9])]:
        datas=[date(p['ano_financeiro'],m,d) for m in meses
               for d in range(1,calendar.monthrange(p['ano_financeiro'],m)[1]+1)]
        uteis=sum(d.weekday()<5 and d not in feriados for d in datas)
        out[est]=dict(meses=meses,dias=len(datas),dias_uteis=uteis,
                      dias_mes=len(datas)/6,dias_uteis_mes=uteis/6,peso_meses=6)
    REGISTRO['calendario_anual']=out
    REGISTRO['convencao_anual']='Dois meses médios sazonais × 6; 365/366 dias, 12 demandas. Agendas representativas, não simulação anual horária. Feriados só quando listados; sem lista, dias úteis=seg–sex.'
    return out


def anual_grupos(contexto,p,curvas,cal,modo='lp',modalidades=('verde',),codigos=('C2','C3')):
    out={}
    for est in ('verao','inverno'):
        q=dict(p,dias_mes=cal[est]['dias_mes'],dias_uteis_mes=cal[est]['dias_uteis_mes'])
        out[est]=grupo(contexto+'/'+est,q,curvas[est],modo=modo,modalidades=modalidades,codigos=codigos)
    return out


def anual_soma(g,mod,cod,chave='opex'):
    rs=[g[e].get((mod,cod),{}) for e in ('verao','inverno')]
    return sum(6*r.get(chave,np.nan) for r in rs) if all(r.get('valido') for r in rs) else np.nan


def classe_capex(c,ref):
    return 'indeterminada' if not _finito(c) else ('viavel' if c>=ref else 'inviavel')


def economia_anual(g,p,mod='verde',a='C2',b='C3'):
    s=anual_soma(g,mod,a)-anual_soma(g,mod,b)
    reg=anual_soma(g,mod,b,'regularizador_rs_mes') if b=='C3' else 0.
    af=fator_anuidade(p['taxa_desconto'],p['vida_bess_anos']);e=p['cap_bess']
    # FOM inclui augmentation: NÃO somar regularizador outra vez à alternativa FOM.
    bruto=s*af/e
    apos_reg=(s-reg)*af/e
    liquido=s*af/(e*(1+p['fom_bess_fracao_ano']*af))
    rs=[r for sem in g.values() for (m,c),r in sem.items() if m==mod and c in (a,b)]
    completo=len(rs)==4 and all(r.get('valido') for r in rs)
    convergido=completo and all(r.get('convergido') for r in rs)
    return dict(economia_anual_R=s,regularizador_anual_R=reg,CAPEX_bruto_RkWh=bruto,
        CAPEX_apos_regularizador_RkWh=apos_reg,CAPEX_liquido_FOM_RkWh=liquido,
        classe_bruta=classe_capex(bruto,p['custo_bess_kwh']),classe_liquida=classe_capex(liquido,p['custo_bess_kwh']),
        banda_liquida=('indeterminada' if not _finito(liquido) else 'viavel_em_toda_faixa' if liquido>=4000 else
                       'inviavel_em_toda_faixa' if liquido<1100 else 'depende_do_CAPEX'),
        aceito=completo,convergido=convergido,
        leitura='estimativa do despacho com regularizador; FOM é alternativa econômica, não soma ao regularizador')


def financeiro_anual(g,p):
    for mod in ('verde','azul'):
        for comp,a,b,vida,cap,custos,unidade in (
            ('FV','C1','C2',p['vida_fv_anos'],p['pot_fv'],p['capex_fv_faixa'],'R$/kWp'),
            ('BESS','C2','C3',p['vida_bess_anos'],p['cap_bess'],p['capex_bess_faixa'],'R$/kWh'),
            ('FV+BESS','C1','C3',p['vida_bess_anos'],1.,[p['pot_fv']*p['custo_pv_kwp']+p['cap_bess']*p['custo_bess_kwh']],'R$ total')):
            s=anual_soma(g,mod,a)-anual_soma(g,mod,b);af=fator_anuidade(p['taxa_desconto'],vida)
            rs=[g[e][mod,c] for e in ('verao','inverno') for c in (a,b)]
            certificado=all(r.get('valido') and r.get('convergido') for r in rs)
            for c in custos:
                investimento=c*cap
                custo_fom=p['fom_bess_fracao_ano']*(investimento if comp=='BESS' else p['cap_bess']*p['custo_bess_kwh'] if comp=='FV+BESS' else 0.)
                for tipo,economia in [('bruto',s),('liquido_FOM_BESS',s-custo_fom)] if comp!='FV' else [('bruto_sem_OM_FV',s)]:
                    FINANCEIRO.append(dict(modalidade=mod,componente=comp,tipo=tipo,
                        convergido=certificado,estados=' / '.join(estado(r) for r in rs),
                        economia_anual_R=economia,CAPEX_unitario=c,unidade=unidade,CAPEX_total_R=investimento,
                        horizonte_anos=vida,payback_anos=investimento/economia if _finito(economia) and economia>0 else np.nan,
                        VPL_R=-investimento+economia*af,CAPEX_bruto_equilibrio=s*af/cap,
                        CAPEX_liquido_equilibrio=s*af/(cap*(1+p['fom_bess_fracao_ano']*af)) if comp=='BESS' else np.nan,
                        CAPEX_apos_regularizador=(s-anual_soma(g,mod,b,'regularizador_rs_mes'))*af/cap if comp=='BESS' else np.nan,
                        observacao='FV+BESS: 10 anos, sem residual FV; O&M FV não estimado' if comp=='FV+BESS' else
                        'FOM inclui reposição; não soma regularizador' if comp=='BESS' else 'bruto, sem O&M FV; dois perfis anuais'))


def fronteira(p,curvas,cal):
    pontos={};confirmados={}
    for h in p['perm_sweep_h']:
        for c in p['credito_sweep']:
            q=_q_perm(p,h);q['multiplicador_credito']=c
            g=anual_grupos(f'fronteira/{h}h/credito{c}/LP',q,curvas,cal)
            v=economia_anual(g,q);pontos[h,c]=(q,v)
            FRONTEIRA.append(dict(permanencia_h=h,credito=c,metodo='LP',confirmado_MILP=False,**v))
            publicar()
    hs=sorted(p['perm_sweep_h']);cs=sorted(p['credito_sweep'])
    selecionados={(hs[0],cs[0]),(hs[0],cs[-1]),(hs[-1],cs[0]),(hs[-1],cs[-1]),(10.,1.)}
    selecionados &= set(pontos)
    def selecionar_vizinhos():
        for i,h in enumerate(hs):
            for j,c in enumerate(cs):
                a=confirmados.get((h,c),pontos[h,c][1])
                viz=[]
                if i+1<len(hs):viz.append((hs[i+1],c))
                if j+1<len(cs):viz.append((h,cs[j+1]))
                for k in viz:
                    b=confirmados.get(k,pontos[k][1])
                    if any(a[cl]!='indeterminada' and b[cl]!='indeterminada' and a[cl]!=b[cl]
                           for cl in ('classe_bruta','classe_liquida')):
                        selecionados.update(((h,c),k))
        # Pontos LP reprovados não são escondidos: tentar confirmação e reportar origem.
        selecionados.update(k for k,(q,v) in pontos.items() if not v['aceito'])
    selecionar_vizinhos()
    while selecionados-set(confirmados):
        for h,c in sorted(selecionados-set(confirmados)):
            q,lp=pontos[h,c]
            g=anual_grupos(f'fronteira/{h}h/credito{c}/MILP',q,curvas,cal,modo='milp')
            mi=economia_anual(g,q);confirmados[h,c]=mi
            row=next(x for x in FRONTEIRA if x['permanencia_h']==h and x['credito']==c)
            row.update(confirmado_MILP=bool(mi['aceito'] and mi['convergido']),
                economia_MILP_anual_R=mi['economia_anual_R'],CAPEX_MILP_bruto_RkWh=mi['CAPEX_bruto_RkWh'],
                CAPEX_MILP_liquido_RkWh=mi['CAPEX_liquido_FOM_RkWh'],
                classe_MILP_bruta=mi['classe_bruta'],classe_MILP_liquida=mi['classe_liquida'],
                diferenca_LP_MILP_anual_R=mi['economia_anual_R']-lp['economia_anual_R'])
            publicar()
        selecionar_vizinhos()
    for tipo in ('bruta','liquida'):
        classes=[(confirmados[k] if k in confirmados and confirmados[k]['aceito'] else v)['classe_'+tipo]
                 for k,(q,v) in pontos.items()]
        completo=all(v['aceito'] and v['convergido'] for q,v in pontos.values()) and all(v['aceito'] and v['convergido'] for v in confirmados.values())
        REGISTRO['H2_'+tipo]=dict(viaveis=classes.count('viavel'),inviaveis=classes.count('inviavel'),
            indeterminadas=classes.count('indeterminada'),
            veredito=('inconclusiva' if not completo or 'indeterminada' in classes else
                      'corroborada nos pontos examinados' if len(set(classes))==2 else 'refutada nos pontos examinados'),
            ressalva='grade aproximada; classes pontuais com gaps reportados; não certifica todo o contínuo entre pontos')


def verificacoes_complementares(p):
    af=fator_anuidade(.10,10)
    dados=[('FIN','economia mensal exigida (R$)',20343.,1500000/(12*af),1.),
           ('FIN','CAPEX* com 24340,49/mês (R$/kWh)',1794.74,24340.49*12*af/1000,.02),
           ('FIN','CAPEX* com 35034,759646/mês (R$/kWh)',2583.28,35034.759646*12*af/1000,.02),
           ('TAR','ultrapassagem D400/P585,73 (R$)',4465.,custo_demanda_verde_fixo(585.73,400,p)[1],1.),
           ('TAR','ultrapassagem D400/P604,28 (R$)',4911.,custo_demanda_verde_fixo(604.28,400,p)[1],1.)]
    # Caso G controlado: 1 kWh CA entregue na ponta, recarregado fora dela, eta_rt=.9025.
    # A extensão de 764 kWh é ORDEM DE GRANDEZA: não teste de igualdade da economia total C2–C3.
    eta2=.95**2;unit=1.81466-.45128/eta2
    m=pulp.LpProblem('Verificacao_G_arbitragem',pulp.LpMinimize)
    ch=pulp.LpVariable('G_carga',lowBound=0,upBound=2)
    dis=pulp.LpVariable('G_descarga',lowBound=0,upBound=1)
    m+=eta2*ch==dis  # fechamento de energia: entrada CA → bateria → saída CA
    m+=.45128*ch+1.81466*(1-dis)
    solver,_=criar_solver(p,p['tempo_limite_s']);ini=time.perf_counter();m.solve(solver)
    certificado=m.status==pulp.LpStatusOptimal and m.sol_status==pulp.LpSolutionOptimal
    obt=1.81466-float(pulp.value(m.objective)) if certificado else np.nan
    REGISTRO_SOLVES.append(dict(id='G',tentativa=1,rotulo='Verificacao G',bloco='verificacao',modo='LP',
        solver='HiGHS',status=pulp.LpStatus.get(m.status),gap_pct=0. if certificado else np.nan,
        convergido=certificado,tempo_s=time.perf_counter()-ini,limite_tempo_s=p['tempo_limite_s'],
        obj=float(pulp.value(m.objective)) if certificado else np.nan,n_var=2,n_bin=0))
    dados.append(('G','economia por 1 kWh CA deslocado (R$)',unit,obt,1e-7))
    for cod,gr,esp,obt,tol in dados:
        VERIFICACOES.append(dict(caso=cod,grandeza=gr,esperado=esp,obtido=obt,erro=obt-esp,
                                 tolerancia=tol,aprovado=abs(obt-esp)<=tol))
    REGISTRO['caso_G_ordem_grandeza']=dict(energia_CA_kWh_dia=764.,dias_ponta=22,
        recarga_rede_fora_Rmes=764*22*unit,
        excedente_FV_credito_zero_Rmes=764*22*1.81466,
        energia_entrada_necessaria_kWh_dia=764/eta2,
        ressalva='Não é previsão de C2-C3: requer excedente/energia suficientes, não inclui demanda, demais dias nem crédito. As perdas incidem na energia de entrada.')


def regularizador(p,curva):
    for h,c in [(10.,1.),(4.,0.),(4.,1.)]:
        q=_q_perm(p,h);q['multiplicador_credito']=c
        a=grupo(f'robustez/{h}h/{c}/padrao',q,curva,codigos=('C2','C3'))
        z=grupo(f'robustez/{h}h/{c}/zero',dict(q,custo_degrad_otim=0.),curva,codigos=('C2','C3'))
        d=diferenca(a['verde','C2'],a['verde','C3']);dz=diferenca(z['verde','C2'],z['verde','C3'])
        reg=a['verde','C3'].get('regularizador_rs_mes',np.nan)
        COMPARACOES.append(dict(bloco=f'robustez/{h}h/{c}',modalidade='verde',comparacao='regularizador',
            economia_Rmes=d,economia_sem_regularizador_Rmes=dz,regularizador_Rmes=reg,
            regularizador_pct=100*reg/d if _finito(d) and d>p['tol_fatura_rs'] else np.nan,
            diferenca_Rmes=dz-d,observacao='robustez pontual; não generalizar a todos os casos'))


def executar_matriz(p,curvas,cal,modo_execucao):
    tip=curvas['inverno'];mods=('verde','azul')
    base=grupo('base',p,tip,modalidades=mods,perfil=True)
    if modo_execucao=='conferencia':
        # Regressões que a rodada anterior revelou: alvo, crédito zero, 6 e 8 h.
        for h in (6.,8.):
            q=_q_perm(p,h);q['multiplicador_credito']=0.
            grupo(f'conferencia/alvo/{h}h',q,tip,modo='lp',codigos=('C2','C3'))
        grupo('conferencia/dois_carregadores',p,tip,N=2,codigos=('C1','C2','C3'))
        grupo('conferencia/P10_observado',p,curvas['p10_observado'],codigos=('C2','C3'),perfil=True)
        anual=anual_grupos('anual',p,curvas,cal,modo='milp',modalidades=mods,codigos=('C1','C2','C3'))
        financeiro_anual(anual,p)
        return
    for f in (20,40,54,60,80,100):grupo(f'frota/{f}',p,tip,f=f,modalidades=mods)
    for n in p['carregadores_lista']:grupo(f'carregadores/{n}',p,tip,N=n)
    for h in p['perm_sweep_h']:grupo(f'permanencia/{h}h',_q_perm(p,h),tip,perfil=h==4)
    for c in p['credito_sweep']:grupo(f'credito/{c}',dict(p,multiplicador_credito=c),tip)
    for pf in ('walz_rudion','noturno_hipotese'):
        for est in ('verao','inverno'):
            grupo(f'perfil/{pf}/{est}',dict(p,perfil_chegada=pf),curvas[est])
    for nome in ('p5_p15','p10_observado'):
        sol=grupo('solar/'+nome,p,curvas[nome],codigos=('C2','C3'),perfil=True)
        COMPARACOES.append(dict(bloco='solar/'+nome,modalidade='verde',comparacao='FV',
            economia_Rmes=diferenca(base['verde','C1'],sol['verde','C2']),
            observacao='C1 do caso-base, mesma frota; sem dependência da geração FV'))
    for nome in ('inverno','p5_p15'):
        for d in p['contrato_sweep_kw']:
            grupo(f'contrato/{nome}/{d}',p,curvas[nome],D=d,modo='lp',codigos=('C2','C3'))
    anual=anual_grupos('anual',p,curvas,cal,modo='milp',modalidades=mods,codigos=('C0','C0b','C1','C2','C3'))
    financeiro_anual(anual,p)
    fronteira(p,curvas,cal)
    for cap in p['capacidades_bess_kwh']:
        for credito in (0.,1.):
            q=_q_perm(p,4.);q.update(cap_bess=cap,multiplicador_credito=credito)
            g=anual_grupos(f'capacidade/{cap}/{credito}',q,curvas,cal,modo='milp')
            DIMENSIONAMENTO.append(dict(capacidade_kWh=cap,potencia_kW=q['pot_bess'],credito=credito,
                ciclos_equivalentes_ano=anual_soma(g,'verde','C3','efc'),
                descarga_BESS_kWh_ano=anual_soma(g,'verde','C3','efc')*cap,
                energia_importada_ponta_kWh_ano=anual_soma(g,'verde','C3','ponta_kwh_mes'),
                energia_util_CA_kWh=cap*(q['soc_max_bess']-q['soc_min_bess'])*q['eta_bess'],**economia_anual(g,q)))
    regularizador(p,tip)
    for seed in range(p['n_realizacoes_monte_carlo']):
        grupo(f'sementes/LP/{seed}',p,tip,seed=seed,modo='lp')
        if seed<3:grupo(f'sementes/MILP/{seed}',p,tip,seed=seed,modo='milp')


def validar(p):
    if not p['alvo_exato']:raise ValueError('Esta versão requer SoC final igual ao alvo.')
    if int(p['n_realizacoes_monte_carlo'])!=p['n_realizacoes_monte_carlo']:raise ValueError('Número de sementes deve ser inteiro.')
    if p['frota_base']!=54 or p['perm_media_h']!=10 or p['cap_bat']!=300:
        raise ValueError('Esta matriz E5 mantém frota-base 54, permanência 10 h e bateria 300 kWh.')
    if p['contrato_modo']!='otimo' or p['regra_ultrapassagem']!='faturavel_max':
        raise ValueError('Mantenha o caso-base redimensionado ao pico e a regra faturavel_max; contratos fixos são sensibilidade.')
    if p['serie_fv_origem']!='NASA_POWER' or p['ocupacao']!='livre':raise ValueError('Fonte/ocupação incompatível com esta metodologia.')
    if not 0<p['eta_carga']<=1 or not 0<p['eta_bess']<=1:raise ValueError('Rendimentos fora de (0,1].')
    if not 0<=p['soc_i_min']<=p['soc_i_max']<=p['soc_f']<=1:raise ValueError('Faixa de SoC inválida.')
    if not 0<=p['soc_min_bess']<p['soc_max_bess']<=1:raise ValueError('Limites do BESS inválidos.')
    for k in ('cap_bess','pot_bess','pot_fv','pot_carregador','n_carregadores','tempo_limite_s','tempo_max_s','limite_total_caso_s'):
        if not _finito(p[k]) or p[k]<=0:raise ValueError('Parâmetro deve ser positivo e finito: '+k)
    if p['tempo_max_s']<p['tempo_limite_s']:raise ValueError('Tempo máximo menor que inicial.')
    if not 0<p['gap_rel']<1 or p['fom_bess_fracao_ano']<0 or p['custo_degrad_otim']<0:raise ValueError('Custos/gap inválidos.')
    if not 0<p['dias_uteis_mes']<=p['dias_mes']:raise ValueError('Calendário mensal inválido.')
    if p['n_realizacoes_monte_carlo']<3:raise ValueError('E5 requer pelo menos três sementes.')
    if p['threads_solver']!=1:raise ValueError('Esta versão fixa uma thread para facilitar reprodutibilidade.')
    if any(x<0 or x>1 for x in p['credito_sweep']):raise ValueError('Crédito fora de [0,1].')
    if p['perm_sweep_h']!=[4.,6.,8.,10.] or p['credito_sweep']!=[0.,.5,1.]:
        raise ValueError('Matriz E5 aprovada: permanência 4/6/8/10 h e crédito 0/0,5/1.')


def imprimir_parametros(p):
    print('\nSIMEA 2027 — V22.3 | parâmetros agrupados')
    for titulo,linhas in [
        ('FROTA',[f"54 veículos/dia; bateria 300 kWh; alvo {p['soc_f']:.0%}; eficiência {p['eta_carga']:.0%}",
                  f"Chegadas: {p['perfil_chegada']}; permanência base 10 h; {p['n_carregadores']} × {p['pot_carregador']:g} kW"]),
        ('FV E BESS',[f"Guarulhos; NASA 2019–2023; FV {p['pot_fv']:g} kWp; PR {p['pr_fv']:.2f}",
                       f"BESS {p['cap_bess']:g} kWh / {p['pot_bess']:g} kW; ciclo fechado, sem uso mínimo imposto"]),
        ('ECONOMIA',["Tarifas históricas E4, Verde e Azul; não orçamento atual", 
                      f"BESS referência {p['custo_bess_kwh']:g} R$/kWh; FOM {100*p['fom_bess_fracao_ano']:g}%/ano (proxy ATB, inclui reposição)",
                      "Regularizador fora da fatura; alternativa FOM não soma regularizador"]),
        ('EXECUÇÃO',[f"HiGHS; gap {100*p['gap_rel']:g}%; teto {p['tempo_max_s']:g} s/tentativa e {p['limite_total_caso_s']:g} s/caso estrito",
                      f"{p['n_realizacoes_monte_carlo']} sementes; saídas consolidadas, sem arquivos por caso"])]:
        print('\n'+titulo)
        for s in linhas:print('  '+s)


def editar(p):
    if input('\nManter os parâmetros padrão? [S/n]: ').strip().lower() in ('','s','sim'):return p
    print('Caso-base físico preservado. Ajustes de execução e custos de sensibilidade:')
    for texto,k,tipo in [('Tempo inicial (s)','tempo_limite_s',float),('Tempo máximo/tentativa (s)','tempo_max_s',float),
                          ('Teto acumulado/caso estrito (s)','limite_total_caso_s',float),('Número de sementes (3 E5, 10 E6)','n_realizacoes_monte_carlo',int),
                          ('FOM anual como fração do CAPEX, inclui reposição','fom_bess_fracao_ano',float)]:
        p[k]=obter_input(texto,p[k],tipo)
    return p


def resumo_final(concluida,p):
    dados=pd.DataFrame(LINHAS)
    counts=dados.estado.value_counts().to_dict() if len(dados) else {}
    REGISTRO['contagem_estados']=counts
    REGISTRO['concluida']=concluida
    REGISTRO['validacao_numerica_pelo_editor']=False
    linhas=['# Leitura da rodada V22.3',
        'EXECUÇÃO CONCLUÍDA' if concluida else 'EXECUÇÃO INTERROMPIDA — arquivos parciais, não usar como rodada final.',
        'Modo: '+REGISTRO['modo']+'. Rodada de conferência não substitui a matriz E5.',
        'Estados (incluem o mesmo caso em blocos distintos): '+str(counts),
        '1. Comece por tabelas/verificacao.csv e auditoria/auditoria.csv.',
        '2. Em tabelas/cenarios.csv, filtre bloco=base: dez linhas nas duas modalidades.',
        '3. VIAVEL_COM_GAP mantém a fatura do incumbente, sem afirmar ótimo. Consulte auditoria/solver.csv.',
        '4. Limites fatura_min são do mínimo da fatura isolada. Em C3 não equivalem ao ótimo de fatura+regularizador.',
        '5. financeiro_anual usa dois semestres representativos, não desempenho anual observado. Baixa geração é estresse.',
        '6. fronteira usa CAPEX anual bruto e líquido; FOM já inclui reposição. Não somar ao regularizador.',
        '7. As figuras têm células N/D e pontos inválidos omitidos explicitamente; não interpolar conclusões através de lacunas.',
        '8. Verificar legibilidade no template antes de publicar. Fonte nominal 8 pt a 17,5 cm.',
        'A execução do editor não foi realizada: esta rodada é evidência produzida pelo usuário.']
    for k in ('H2_bruta','H2_liquida'):
        if k in REGISTRO:linhas.append(k+': '+str(REGISTRO[k]))
    caminho('LEIA_PRIMEIRO.md').write_text('\n\n'.join(linhas)+'\n',encoding='utf8')


def main():
    parser=argparse.ArgumentParser(description='SIMEA 2027 V22.3 — execução pelo autor, saídas consolidadas')
    parser.add_argument('--modo',choices=['verificacao','conferencia','e5'],help='sem argumento: menu interativo')
    parser.add_argument('--config',type=Path,help='JSON parcial de parâmetros; alterações registradas')
    parser.add_argument('--dados',type=Path,help='pasta do cache NASA com CSV e JSON originais')
    parser.add_argument('--saida',type=Path,default=Path(__file__).resolve().parent/'resultados')
    args=parser.parse_args();p=parametros_padrao()
    if args.config:
        alter=json.loads(args.config.read_text(encoding='utf8-sig'))
        if not isinstance(alter,dict) or set(alter)-set(p):raise ValueError('Configuração contém chaves desconhecidas.')
        p.update(alter)
    p['arquivo_serie_fv']=str((args.dados or Path(__file__).resolve().parent/'dados_nasa_v22').resolve())
    imprimir_parametros(p)
    modo=args.modo
    if modo is None:
        editar(p)
        modo=input('\nRodada [verificacao/conferencia/e5; padrão conferencia]: ').strip().lower() or 'conferencia'
    if modo not in ('verificacao','conferencia','e5'):raise ValueError('Modo inválido.')
    validar(p)
    dest=args.saida.resolve()/('E5_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')+'_'+modo)
    dest.mkdir(parents=True,exist_ok=False);_PASTA['saida']=str(dest)
    shutil.copyfile(__file__,caminho('auditoria/codigo_executado.py'))
    REGISTRO.update(versao=VERSAO,codigo_sha256=_hash_codigo(),origem_sha256=CODIGO_ORIGEM_SHA256,
        ambiente=versoes_ambiente(),parametros=copy.deepcopy(p),modo=modo,inicio=datetime.now().astimezone().isoformat(),
        concluida=False,fonte_FOM='https://atb.nrel.gov/electricity/2024b/commercial_battery_storage',
        nota_FOM='2,5% CAPEX/ano inclui augmentation; proxy internacional com horizonte original 15 anos, aqui 10; sensibilidade, não cotação local')
    print('\nSaída:',dest)
    inicio=time.perf_counter();concluida=False
    try:
        if not executar_verificacao(p):raise RuntimeError('A–F reprovados: interromper e revisar verificacao.csv.')
        verificacoes_complementares(p)
        if not all(x['aprovado'] for x in VERIFICACOES):raise RuntimeError('Verificação complementar reprovada.')
        if modo!='verificacao':
            curvas=preparar_entradas(p);cal=calendario_anual(p)
            executar_matriz(p,curvas,cal,modo)
            gerar_figuras(p)
        concluida=True
    except BaseException as exc:
        REGISTRO['erro']=dict(tipo=type(exc).__name__,mensagem=str(exc),rastro=traceback.format_exc())
        print('\nRodada interrompida. Consulte auditoria/registro_execucao.json.')
        raise
    finally:
        REGISTRO.update(fim=datetime.now().astimezone().isoformat(),tempo_total_s=time.perf_counter()-inicio)
        resumo_final(concluida,p);publicar()
    print('\nRodada concluída. Leia LEIA_PRIMEIRO.md; conclusão não significa aprovação de todos os casos.')


# -----------------------------------------------------------------------------
# Figuras para leitura: pontos com gap aparecem vazados; sem pontos inválidos em zero.
# -----------------------------------------------------------------------------
CORES={'C0':'#555555','C0b':'#999999','C1':'#0072B2','C2':'#d89000','C3':'#009E73'}


def figura(nome,n=2):
    plt.rcParams.update({'font.size':8,'axes.labelsize':8,'legend.fontsize':8,'axes.titlesize':9})
    fig,axs=plt.subplots(1,n,figsize=(17.5/2.54,8.5/2.54),squeeze=False,layout='constrained')
    return fig,axs.ravel(),nome


def guardar(fig,nome):
    fig.savefig(caminho('figuras/'+nome+'.png'),dpi=300)
    plt.close(fig)
    REGISTRO.setdefault('figuras',[]).append(nome+'.png')


def curvas_cenarios(ax,df,x,y):
    for cod in CORES:
        g=df[(df.cenario==cod)&(df.modalidade=='verde')].sort_values(x)
        if g.empty:continue
        ax.plot(g[x],g[y],'-',color=CORES[cod],linewidth=1,label=cod)
        ok=g.estado.ne('VIAVEL_COM_GAP')
        ax.scatter(g.loc[ok,x],g.loc[ok,y],s=12,color=CORES[cod])
        ax.scatter(g.loc[~ok,x],g.loc[~ok,y],s=18,facecolors='none',edgecolors=CORES[cod])
    ax.grid(alpha=.2)


def gerar_figuras(p):
    d=pd.DataFrame(LINHAS);cmp=pd.DataFrame(COMPARACOES);pr=pd.DataFrame(PERFIS)
    fig,axs,nome=figura('01_operacao_base')
    for cod in CORES:
        g=pr[(pr.bloco=='base')&(pr.modalidade=='verde')&(pr.cenario==cod)]
        if len(g):axs[0].step(g.hora,g.rede_kW,where='post',color=CORES[cod],label=cod,
                             linestyle='--' if cod in ('C0b','C2') else '-')
    axs[0].axvspan(18,21,alpha=.1,color='gray');axs[0].set(xlabel='Hora LST',ylabel='Importação (kW)')
    axs[0].legend(ncol=3,loc='upper left')
    g=pr[(pr.bloco=='base')&(pr.modalidade=='verde')&(pr.cenario=='C3')]
    for col,label in [('FV_kW','FV'),('BESS_carga_kW','Carga BESS'),('BESS_descarga_kW','Descarga BESS')]:
        if len(g):axs[1].step(g.hora,g[col],where='post',label=label)
    axs[1].set(xlabel='Hora LST',ylabel='Potência (kW)');axs[1].legend(loc='upper right');guardar(fig,nome)
    for bloco,nome,x,xlabel,ys in [
        ('frota/','02_frota','frota','Veículos/dia',[('pico_kW','Pico (kW)'),('fatura_Rmes','Fatura (R$/mês)')]),
        ('carregadores/','03_carregadores','carregadores','Carregadores',[('deficit_kWh_dia','Déficit (kWh/dia)'),('fatura_Rmes','Fatura (R$/mês)')]),
        ('permanencia/','04_permanencia','permanencia_nominal_h','Permanência nominal (h)',[('energia_ponta_kWh_dia','Energia na ponta (kWh/dia)'),('fatura_Rmes','Fatura (R$/mês)')])]:
        g=d[d.bloco.str.startswith(bloco)]
        if g.empty:continue
        fig,axs,n=figura(nome)
        for ax,(y,yl) in zip(axs,ys):curvas_cenarios(ax,g,x,y);ax.set(xlabel=xlabel,ylabel=yl)
        axs[0].legend(ncol=3,loc='upper right')
        fig.suptitle('Pontos vazados: solução viável com gap. Valores reprovados não compõem faturas.',fontsize=8)
        guardar(fig,n)
    if FRONTEIRA:
        grid=pd.DataFrame(FRONTEIRA);hs=sorted(grid.permanencia_h.unique());cs=sorted(grid.credito.unique())
        fig,axs,n=figura('05_fronteira_anual')
        for ax,col,title in zip(axs,['CAPEX_bruto_RkWh','CAPEX_liquido_FOM_RkWh'],['Bruto','Líquido de FOM/reposição']):
            z=grid.pivot(index='credito',columns='permanencia_h',values=col).reindex(index=cs,columns=hs).to_numpy()
            cm=plt.get_cmap('cividis').copy();cm.set_bad('#dddddd')
            im=ax.imshow(np.ma.masked_invalid(z),origin='lower',aspect='auto',cmap=cm,vmin=0)
            finite=z[np.isfinite(z)];mx=max(float(finite.max()),1.) if finite.size else 1.
            for i,c in enumerate(cs):
                for j,h in enumerate(hs):
                    val=z[i,j];ax.text(j,i,f'{val:.0f}' if np.isfinite(val) else 'N/D',ha='center',va='center',
                                       color='black' if not np.isfinite(val) or val>mx*.55 else 'white')
                    row=grid[(grid.permanencia_h==h)&(grid.credito==c)].iloc[0]
                    if row.confirmado_MILP:ax.plot(j+.32,i+.30,'o',ms=3,color='red')
            if finite.size:
                levels=[a for a in [1100,1500,4000] if finite.min()<a<finite.max()]
                if levels:
                    cont=ax.contour(np.ma.masked_invalid(z),levels=levels,colors='white',linewidths=.8)
                    ax.clabel(cont,fmt='%d',fontsize=7)
            ax.set_xticks(range(len(hs)),hs);ax.set_yticks(range(len(cs)),cs)
            ax.set(xlabel='Permanência (h)',ylabel='Multiplicador do crédito',title=title)
            fig.colorbar(im,ax=ax,label='CAPEX* (R$/kWh)',shrink=.8)
        fig.suptitle('Grade LP anual; ponto vermelho: conferido em MILP. Contornos aproximados da faixa de CAPEX.',fontsize=8)
        guardar(fig,n)
    g=d[d.bloco.str.startswith('contrato/')]
    if len(g):
        fig,axs,n=figura('06_contrato')
        for serie,cor in [('inverno','#0072B2'),('p5_p15','#D55E00')]:
            for cod,ls in [('C2','--'),('C3','-')]:
                z=g[g.bloco.str.startswith('contrato/'+serie+'/')&g.cenario.eq(cod)].sort_values('contrato_kW')
                axs[0].plot(z.contrato_kW,z.fatura_Rmes,ls+'o',ms=3,color=cor,label=serie+' '+cod)
            z=cmp[cmp.bloco.str.startswith('contrato/'+serie+'/')&cmp.comparacao.eq('BESS')].copy()
            z['D']=z.bloco.str.rsplit('/',n=1).str[-1].astype(float);z=z.sort_values('D')
            axs[1].plot(z.D,z.economia_Rmes,'o-',ms=3,color=cor,label=serie)
        axs[0].set(xlabel='Contrato (kW)',ylabel='Fatura (R$/mês)');axs[1].set(xlabel='Contrato (kW)',ylabel='Economia BESS (R$/mês)')
        for ax in axs:ax.legend(fontsize=7);ax.grid(alpha=.2)
        fig.suptitle('Ocupação relaxada; demanda faturável=max(contrato, medida).',fontsize=8);guardar(fig,n)
    g=d[d.bloco.str.startswith('sementes/LP/')]
    if len(g):
        fig,axs,n=figura('07_sementes')
        for cod in ('C0','C1','C3'):
            z=g[g.cenario.eq(cod)];axs[0].scatter(z.seed,z.pico_kW,color=CORES[cod],label=cod,s=15)
        for met,marker in [('LP','o'),('MILP','x')]:
            z=cmp[cmp.bloco.str.startswith('sementes/'+met+'/')&cmp.comparacao.eq('BESS')]
            xx=z.bloco.str.rsplit('/',n=1).str[-1].astype(int)
            axs[1].scatter(xx,z.economia_Rmes,marker=marker,label=met,s=22)
        axs[0].set(xlabel='Semente',ylabel='Pico (kW)');axs[1].set(xlabel='Semente',ylabel='Economia BESS (R$/mês)')
        for ax in axs:ax.legend()
        guardar(fig,n)
    if (d.bloco=='solar/p10_observado').any():
        fig,axs,n=figura('08_baixa_geracao')
        x=np.arange(3);blocos=['base','solar/p5_p15','solar/p10_observado']
        for cod,off in [('C2',-.18),('C3',.18)]:
            vals=[d[(d.bloco==b)&(d.modalidade=='verde')&(d.cenario==cod)].pico_kW.iloc[0] for b in blocos]
            axs[0].bar(x+off,vals,.36,label=cod,color=CORES[cod])
        for tipo,off,cor in [('FV',-.18,CORES['C2']),('BESS',.18,CORES['C3'])]:
            vals=[]
            for b in blocos:
                z=cmp[(cmp.bloco==b)&(cmp.modalidade=='verde')&(cmp.comparacao==tipo)]
                vals.append(z.economia_Rmes.iloc[0] if len(z) else np.nan)
            axs[1].bar(x+off,vals,.36,label=tipo,color=cor)
        for ax in axs:ax.set_xticks(x,['Típico','Média\nP5–P15','Dia\nP10']);ax.legend(loc='lower center',bbox_to_anchor=(.5,1),ncol=2)
        axs[0].set_ylabel('Pico (kW)');axs[1].set_ylabel('Economia (R$/mês)');guardar(fig,n)
    if DIMENSIONAMENTO:
        fig,axs,n=figura('09_capacidade_BESS');dd=pd.DataFrame(DIMENSIONAMENTO)
        for c,g in dd.groupby('credito'):
            g=g.sort_values('capacidade_kWh')
            axs[0].plot(g.capacidade_kWh,g.economia_anual_R,'o-',label=f'Crédito {c:g}')
            axs[1].plot(g.capacidade_kWh,g.CAPEX_liquido_FOM_RkWh,'o-',label=f'Crédito {c:g}')
        axs[0].set(xlabel='Capacidade nominal (kWh)',ylabel='Economia bruta anual (R$/ano)')
        axs[1].set(xlabel='Capacidade nominal (kWh)',ylabel='CAPEX* líquido (R$/kWh)')
        axs[1].axhspan(1100,4000,color='gray',alpha=.15,label='Faixa 1.100–4.000')
        axs[1].axhline(1500,color='black',ls='--',lw=1)
        for ax in axs:ax.legend(fontsize=7);ax.grid(alpha=.2)
        fig.suptitle('Permanência 4 h; potência BESS fixa em 500 kW.',fontsize=8);guardar(fig,n)


if __name__=='__main__':
    main()

