from cobra import Model, Reaction, Metabolite
import pandas as pd

# 1. Khởi tạo mô hình
model = Model("ecoli_acetate_model_final")
model.compartments['c'] = 'cytosol'
model.compartments['e'] = 'extracellular'

# 2. Định nghĩa Metabolites
# ---- Nội bào ----
_13dpg = Metabolite("13dpg_c", compartment='c')
_3pg = Metabolite("3pg_c", compartment='c')
_2pg = Metabolite("2pg_c", compartment='c')
g3p = Metabolite("g3p_c", compartment='c')
dhap  = Metabolite("dhap_c", compartment='c')
fdp = Metabolite("fdp_c", compartment='c')
g6p   = Metabolite("g6p_c", compartment='c')
f6p   = Metabolite("f6p_c", compartment='c')
pep   = Metabolite("pep_c", compartment='c')
pyr   = Metabolite("pyr_c", compartment='c')
ru5p = Metabolite("ru5p__D_c", compartment='c')
xu5p = Metabolite("xu5p__D_c", compartment='c')
r5p  = Metabolite("r5p_c", compartment='c')
e4p  = Metabolite("e4p_c", compartment='c')
s7p  = Metabolite("s7p_c", compartment='c')
icit = Metabolite("icit_c", compartment='c')
akg  = Metabolite("akg_c", compartment='c')
succ = Metabolite("succ_c", compartment='c')
fum  = Metabolite("fum_c", compartment='c')
mal  = Metabolite("mal__L_c", compartment='c')
oaa = Metabolite("oaa_c", compartment='c')
cit = Metabolite("cit_c", compartment='c')
accoa = Metabolite("accoa_c", compartment='c')
glx = Metabolite("glx_c", compartment='c')
ac = Metabolite("ac_c", compartment='c')
ser_L   = Metabolite("ser__L_c", compartment='c')
gly     = Metabolite("gly_c", compartment='c')
glu_L = Metabolite("glu__L_c", compartment='c')
thr_L = Metabolite("thr__L_c", compartment='c')
ile_L = Metabolite("ile__L_c", compartment='c')
val_L = Metabolite("val__L_c", compartment='c')
leu_L = Metabolite("leu__L_c", compartment='c')
phe_L = Metabolite("phe__L_c", compartment='c')
nad   = Metabolite("nad_c", compartment='c')
nadp   = Metabolite("nadp_c", compartment='c')
nadh  = Metabolite("nadh_c", compartment='c')
atp   = Metabolite("atp_c", compartment='c')
adp   = Metabolite("adp_c", compartment='c')
pi    = Metabolite("pi_c", compartment='c')
h    = Metabolite("h_c", compartment='c')
h2o  = Metabolite("h2o_c", compartment='c')
nadph = Metabolite("nadph_c", compartment='c')
co2   = Metabolite("co2_c", compartment='c')
nh4 = Metabolite("nh4_c", compartment='c')

# ---- Ngoại bào ----
ac_e = Metabolite("ac_e", compartment = "e")
co2_e = Metabolite("co2_e", compartment = "e")
nh4_e = Metabolite("nh4_e", compartment="e")
pi_e = Metabolite("pi_e", compartment="e")
h_e = Metabolite("h_e", compartment="e")
h2o_e = Metabolite("h2o_e", compartment="e") 

# 3. Định nghĩa các phản ứng trung tâm
reactions_to_add = []
reaction_specs = [
    ("ATPS", {nadh: -1, adp: -2, pi: -2, h:-2, nad: 1, atp: 2, h2o: 1}, 0, 1000),
    ("ACS", {ac: -1, atp: -1, accoa: 1, adp: 1}, 0, 1000),
    ("CS", {oaa: -1, accoa: -1, h2o: -1, cit: 1, h: 1}, 0, 1000),
    ("ACONT", {cit: -1, icit: 1}, -1000, 1000),
    ("ICDH", {icit: -1, nad: -1, akg: 1, co2: 1, nadh: 1}, 0, 1000),
    ("ICL", {icit: -1, glx: 1, succ: 1}, 0, 1000),
    ("AKGDH", {akg: -1, nad: -1, succ: 1, co2: 1, nadh: 1}, 0, 1000),
    ("SUCD", {succ: -1, nad: -1, fum: 1, nadh: 1}, -1000, 1000),
    ("MALS", {glx: -1, accoa: -1, h2o: -1, mal: 1, h: 1}, 0, 1000),
    ("FUM", {fum: -1, h2o: -1, mal: 1}, -1000, 1000),
    ("MDH", {mal: -1, nad: -1, oaa: 1, nadh: 1}, -1000, 1000),
    ("ME", {mal: -1, nadp: -1, pyr: 1, co2: 1, nadph: 1}, 0, 1000),
    ("PDH", {pyr: -1, nad: -1, accoa: 1, co2: 1, nadh: 1}, 0, 1000),
    ("PCK", {oaa: -1, atp: -1, pep: 1, co2: 1, adp: 1}, -1000, 1000),
    ("ENO", {_2pg: 1, pep: -1, h2o: -1}, -1000, 1000),
    ("PGM", {_3pg: 1, _2pg: -1}, -1000, 1000),
    ("PGK", {_13dpg: 1, adp: 1, _3pg: -1, atp: -1}, -1000, 1000),
    ("GAPD", {g3p: 1, pi: 1, nad: 1, _13dpg: -1, nadh: -1, h: -1}, -1000, 1000),
    ("TPI", {dhap: -1, g3p: 1}, -1000, 1000),
    ("FBA", {fdp: 1, g3p: -1, dhap: -1}, -1000, 1000),
    ("FBP", {fdp: -1, h2o: -1, f6p: 1, pi: 1}, 0, 1000),
    ("PGI", {f6p: -1, g6p: 1}, -1000, 1000),
    ("GND", {g6p: -1, nadp: -2, h2o: -1, ru5p: 1, co2: 1, nadph: 2, h: 2}, 0, 1000),
    ("RPI", {ru5p: -1, r5p: 1}, -1000, 1000),
    ("RPE", {ru5p: -1, xu5p: 1}, -1000, 1000),
    ("TKT1", {r5p: -1, xu5p: -1, g3p: 1, s7p: 1}, -1000, 1000),
    ("TALA", {g3p: -1, s7p: -1, e4p: 1, f6p: 1}, -1000, 1000),
    ("TKT2", {xu5p: -1, e4p: -1, f6p: 1, g3p: 1}, -1000, 1000),
    ("SER_synth", {_3pg: -1, nad: -1, glu_L: -1, ser_L: 1, nadh: 1, akg: 1, h: 1}, 0, 1000),
    ("GLY_synth", {ser_L: -1, gly: 1, co2: 1, nh4: 1}, 0, 1000),
    ("OAA_to_THR", {oaa: -1, thr_L: 1}, 0, 1000),
    ("AKG_to_GLU", {akg: -1, nadph: -1, nh4: -1, h: -1, glu_L: 1, nadp: 1, h2o: 1}, 0, 1000),
    ("OAA_PYR_to_ILE", {oaa: -1, pyr: -1, ile_L: 1}, 0, 1000),
    ("2PYR_to_VAL", {pyr: -2, val_L: 1}, 0, 1000),
    ("AcCoA_2PYR_to_LEU", {accoa: -1, pyr: -2, leu_L: 1}, 0, 1000),
    ("E4P_2PEP_to_PHE", {e4p: -1, pep: -2, phe_L: 1}, 0, 1000)
]
for r_id, stoich, lb, ub in reaction_specs:
    reaction = Reaction(r_id)
    reaction.add_metabolites(stoich)
    reaction.lower_bound = lb
    reaction.upper_bound = ub
    reactions_to_add.append(reaction)
model.add_reactions(reactions_to_add)

# 4. Phản ứng Sinh khối
biomass_precursors = {glu_L: -1.0, thr_L: -1.0, ile_L: -1.0, val_L: -1.0, leu_L: -1.0, phe_L: -1.0, ser_L: -1.0, gly: -1.0, r5p: -1.0, e4p: -1.0, atp: -10}
biomass_reaction = Reaction("BIOMASS_Ecoli_acetate_simple", name="Simplified biomass", lower_bound=0, upper_bound=1000)
biomass_reaction.add_metabolites(biomass_precursors)
model.add_reactions([biomass_reaction])

# 5. Phản ứng Biên
boundary_reactions = []
ex_ac = Reaction("EX_ac_e", lower_bound=-100.0); ex_ac.add_metabolites({ac_e: -1})
ac_tex = Reaction("ACtex", lower_bound=-1000.0, upper_bound=1000.0); ac_tex.add_metabolites({ac_e: -1.0, ac: 1.0})
ex_nh4 = Reaction("EX_nh4_e", lower_bound=-1000.0); ex_nh4.add_metabolites({nh4_e: -1})
nh4_tex = Reaction("NH4tex", lower_bound=-1000.0, upper_bound=1000.0); nh4_tex.add_metabolites({nh4_e: -1.0, nh4: 1.0})
ex_pi = Reaction("EX_pi_e", lower_bound=-1000.0); ex_pi.add_metabolites({pi_e: -1})
pi_tex = Reaction("PItex", lower_bound=-1000.0, upper_bound=1000.0); pi_tex.add_metabolites({pi_e: -1.0, pi: 1.0})
ex_h = Reaction("EX_h_e", lower_bound=-1000.0, upper_bound=1000.0); ex_h.add_metabolites({h_e: -1})
h_tex = Reaction("Htex", lower_bound=-1000.0, upper_bound=1000.0); h_tex.add_metabolites({h_e: -1.0, h: 1.0})
ex_h2o = Reaction("EX_h2o_e", lower_bound=-1000.0, upper_bound=1000.0); ex_h2o.add_metabolites({h2o_e: -1})
h2o_tex = Reaction("H2Otex", lower_bound=-1000.0, upper_bound=1000.0); h2o_tex.add_metabolites({h2o_e: -1.0, h2o: 1.0})
ex_co2 = Reaction("EX_co2_e", lower_bound=0.0, upper_bound=1000.0); ex_co2.add_metabolites({co2_e: -1})
co2_tex = Reaction("CO2tex", lower_bound=-1000.0, upper_bound=1000.0); co2_tex.add_metabolites({co2: -1.0, co2_e: 1.0})
atpm = Reaction("ATPM", name="ATP maintenance", lower_bound=1.0, upper_bound=1000.0); atpm.add_metabolites({atp: -1, h2o: -1, adp: 1, pi: 1, h: 1})
boundary_reactions.extend([ex_ac, ac_tex, ex_nh4, nh4_tex, ex_pi, pi_tex, ex_h, h_tex, ex_h2o, h2o_tex, ex_co2, co2_tex, atpm])
model.add_reactions(boundary_reactions)

# 6. Đặt Hàm mục tiêu và Chạy Mô phỏng
model.objective = "BIOMASS_Ecoli_acetate_simple"
solution = model.optimize()

# 7. Phân tích kết quả
print(f"Trạng thái: {solution.status}")
if solution.status == 'optimal':
    print(f"Giá trị mục tiêu (Tốc độ tăng trưởng tương đối): {solution.objective_value:.4f}")
    print("\n--- Báo cáo Thông lượng (Flux Report) ---")
    print(model.summary()) 
else:
    print("Mô hình không tìm thấy lời giải tối ưu.")