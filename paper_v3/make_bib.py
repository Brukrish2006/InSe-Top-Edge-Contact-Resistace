"""Insert the bibliography into InSe_v3.tex in order of first citation."""
import re, sys
B = {
"refceferino": r"A.~Ceferino, K.W.~Song, S.J.~Magorrian, V.~Z\'olyomi and V.I.~Fal'ko, \textit{Crossover from weakly indirect to direct excitons in atomically thin films of InSe}, Phys.\ Rev.\ B \textbf{101}, 245432 (2020), \doi{10.1103/PhysRevB.101.245432}.",
"refpasquale23": r"G.~Pasquale \textit{et al.}, \textit{Electrical detection of the flat-band dispersion in van der Waals field-effect structures}, Nat.\ Nanotechnol.\ \textbf{18}, 1416 (2023), \doi{10.1038/s41565-023-01489-x}.",
"refsuryavanshi": r"S.V.~Suryavanshi, B.~Magyari-K\"ope, P.~Lim, C.~McClellan, K.K.H.~Smithe, C.D.~English and E.~Pop, \textit{Improving electric contacts to two-dimensional semiconductors}, arXiv:2105.10792 (2021).",
"refpal": r"A.~Pal \textit{et al.}, \textit{Characterization and closed-form modeling of edge/top/hybrid metal--2D semiconductor contacts}, IEDM Tech.\ Dig., 28.5.1 (2022), \doi{10.1109/IEDM45625.2022.10019374}.",
"refdeylgat22": r"E.~Deylgat, E.~Chen, M.V.~Fischetti, B.~Sor\'ee and W.G.~Vandenberghe, \textit{Image-force barrier lowering in top- and side-contacted two-dimensional materials}, Solid-State Electron.\ \textbf{198}, 108458 (2022), \doi{10.1016/j.sse.2022.108458}.",
"refdeylgat23": r"E.~Deylgat, E.~Chen, B.~Sor\'ee and W.G.~Vandenberghe, \textit{Quantum transport study of contact resistance of edge- and top-contacted two-dimensional materials}, Proc.\ SISPAD (2023), \doi{10.23919/SISPAD57422.2023.10319537}.",
"refdeylgat26": r"E.~Deylgat, E.~Chen, B.~Sor\'ee and W.G.~Vandenberghe, \textit{Comparing contact resistance of edge-, top-, and hybrid-contacted two-dimensional materials}, Phys.\ Rev.\ Appl.\ \textbf{25}, 044057 (2026), \doi{10.1103/ddpx-1tt4}.",
"refchen25": r"Z.~Chen \textit{et al.}, \textit{Hybrid-contact engineering for enhanced p-type performance in WSe$_2$ field-effect transistors}, Nanotechnology \textbf{36} (2025), \doi{10.1088/1361-6528/ade15c}.",
"refseok": r"Y.~Seok, S.~Moon, M.~Kim, D.~Yeom, J.~Hwang, K.~Watanabe, T.~Taniguchi and K.~Lee, \textit{Cryogenic hole transport in 2D indium selenide}, ACS Nano \textbf{20}, 24439 (2026), \doi{10.1021/acsnano.6c07680}.",
"refgonzalez": r"S.~Gonzalez-Munoz, K.~Agarwal, E.G.~Castanon, Z.R.~Kudrynskyi, Z.D.~Kovalyuk, J.~Spi\`ece, O.~Kazakova, A.~Patan\`e and O.V.~Kolosov, \textit{Direct measurements of anisotropic thermal transport in $\gamma$-InSe nanolayers via cross-sectional scanning thermal microscopy}, Adv.\ Mater.\ Interfaces \textbf{10}, 2300081 (2023), \doi{10.1002/admi.202300081}.",
"reftung": r"R.T.~Tung, \textit{The physics and chemistry of the Schottky barrier height}, Appl.\ Phys.\ Rev.\ \textbf{1}, 011304 (2014), \doi{10.1063/1.4858400}.",
"refkantner": r"M.~Kantner, \textit{Generalized Scharfetter--Gummel schemes for electro-thermal transport in degenerate semiconductors using the Kelvin formula for the Seebeck coefficient}, J.\ Comput.\ Phys.\ \textbf{402}, 109091 (2020), \doi{10.1016/j.jcp.2019.109091}.",
"refdevsim": r"J.~S\'anchez, \textit{DEVSIM: a TCAD semiconductor device simulator}, J.\ Open Source Softw.\ \textbf{7}, 3898 (2022), \doi{10.21105/joss.03898}.",
"refsg": r"D.L.~Scharfetter and H.K.~Gummel, \textit{Large-signal analysis of a silicon Read diode oscillator}, IEEE Trans.\ Electron Devices \textbf{16}, 64 (1969), \doi{10.1109/T-ED.1969.16566}.",
"refnilsson": r"N.G.~Nilsson, \textit{Empirical approximations for the Fermi energy in a semiconductor with parabolic bands}, Appl.\ Phys.\ Lett.\ \textbf{33}, 653 (1978), \doi{10.1063/1.90452}.",
"refabdel": r"D.~Abdel, P.~Farrell and J.~Fuhrmann, \textit{Assessing the quality of the excess chemical potential flux scheme for degenerate semiconductor device simulation}, Opt.\ Quantum Electron.\ \textbf{53} (2021), \doi{10.1007/s11082-021-02803-4}.",
"reffarrell": r"P.~Farrell, J.~Moatti, M.~O'Donovan, S.~Schulz and T.~Koprucki, \textit{Importance of satisfying thermodynamic consistency in optoelectronic device simulations for high carrier densities}, Opt.\ Quantum Electron.\ \textbf{55}, 978 (2023), \doi{10.1007/s11082-023-05234-5}.",
"refpasquale": r"G.~Pasquale \textit{et al.}, \textit{Flat-band-induced many-body interactions and exciton complexes in a layered semiconductor}, Nano Lett.\ \textbf{22}, 8883 (2022), \doi{10.1021/acs.nanolett.2c02965}.",
"refhamer": r"M.J.~Hamer \textit{et al.}, \textit{Indirect to direct gap crossover in two-dimensional InSe revealed by angle-resolved photoemission spectroscopy}, ACS Nano \textbf{13}, 2136 (2019), \doi{10.1021/acsnano.8b08726}.",
"refhu": r"S.~Hu \textit{et al.}, \textit{Reconfigurable InSe electronics with van der Waals integration}, Adv.\ Electron.\ Mater.\ \textbf{8}, 2101176 (2022), \doi{10.1002/aelm.202101176}.",
"refwei": r"C.~Wei \textit{et al.}, \textit{h-BN-assisted metal contact transfer to InSe for two-dimensional multifunctional electronic devices}, Chin.\ Phys.\ Lett.\ \textbf{41}, 128501 (2024), \doi{10.1088/0256-307X/41/12/128501}.",
"refband": r"D.A.~Bandurin \textit{et al.}, \textit{High electron mobility, quantum Hall effect and anomalous optical response in atomically thin InSe}, Nat.\ Nanotechnol.\ \textbf{12}, 223 (2017), \doi{10.1038/nnano.2016.242}.",
"refsuch": r"S.~Sucharitakul, N.J.~Goble, U.R.~Kumar, R.~Sankar, Z.A.~Bogorad, F.-C.~Chou, Y.-T.~Chen and X.P.A.~Gao, \textit{Intrinsic electron mobility exceeding $10^3$~cm$^2$/(V\,s) in multilayer InSe FETs}, Nano Lett.\ \textbf{15}, 3815 (2015), \doi{10.1021/acs.nanolett.5b00493}.",
"reffeng14": r"W.~Feng, W.~Zheng, W.~Cao and P.~Hu, \textit{Back gated multilayer InSe transistors with enhanced carrier mobilities via the suppression of carrier scattering from a dielectric interface}, Adv.\ Mater.\ \textbf{26}, 6587 (2014), \doi{10.1002/adma.201402427}.",
"refjiang": r"J.~Jiang, L.~Xu, C.~Qiu and L.-M.~Peng, Nature \textbf{616}, 470 (2023), \doi{10.1038/s41586-023-05819-w}.",
"refqin": r"B.~Qin \textit{et al.}, Science \textbf{389}, 299 (2025), \doi{10.1126/science.adu3803}.",
"refho": r"P.-H.~Ho \textit{et al.}, \textit{High-mobility InSe transistors: the role of surface oxides}, ACS Nano \textbf{11}, 7362 (2017), \doi{10.1021/acsnano.7b03531}.",
"refhuang": r"Y.-T.~Huang \textit{et al.}, \textit{High-performance InSe transistors with ohmic contact enabled by nonrectifying barrier-type indium electrodes}, ACS Appl.\ Mater.\ Interfaces \textbf{10}, 33450 (2018), \doi{10.1021/acsami.8b10576}.",
"reflin": r"C.-Y.~Lin, R.K.~Ulaganathan, R.~Sankar and F.-C.~Chou, \textit{Ambipolar field-effect transistors by few-layer InSe with asymmetry contact metals}, AIP Adv.\ \textbf{7}, 075314 (2017), \doi{10.1063/1.4995589}.",
"refkim24": r"M.~Kim \textit{et al.}, \textit{Superior p-type switching in InSe nanosheets for complementary multifunctional systems}, Nano Lett.\ \textbf{24}, 16090 (2024), \doi{10.1021/acs.nanolett.4c04624}.",
"refsong": r"S.~Song \textit{et al.}, \textit{Indium selenides for next-generation low-power computing devices}, Nat.\ Rev.\ Electr.\ Eng.\ \textbf{3}, 185 (2026), \doi{10.1038/s44287-025-00251-w}.",
"refallain": r"A.~Allain, J.~Kang, K.~Banerjee and A.~Kis, \textit{Electrical contacts to two-dimensional semiconductors}, Nat.\ Mater.\ \textbf{14}, 1195 (2015), \doi{10.1038/nmat4452}.",
"refschulman": r"D.S.~Schulman, A.J.~Arnold and S.~Das, \textit{Contact engineering for 2D materials and devices}, Chem.\ Soc.\ Rev.\ \textbf{47}, 3037 (2018), \doi{10.1039/C7CS00828G}.",
"refwangchh": r"Y.~Wang and M.~Chhowalla, \textit{Making clean electrical contacts on 2D transition metal dichalcogenides}, Nat.\ Rev.\ Phys.\ \textbf{4}, 101 (2022), \doi{10.1038/s42254-021-00389-0}.",
"refkim17": r"C.~Kim \textit{et al.}, \textit{Fermi level pinning at electrical metal contacts of monolayer molybdenum dichalcogenides}, ACS Nano \textbf{11}, 1588 (2017), \doi{10.1021/acsnano.6b07159}.",
"refkang": r"J.~Kang, W.~Liu, D.~Sarkar, D.~Jena and K.~Banerjee, Phys.\ Rev.\ X \textbf{4}, 031005 (2014), \doi{10.1103/PhysRevX.4.031005}.",
"refliusa": r"Y.~Liu, P.~Stradins and S.-H.~Wei, \textit{Van der Waals metal-semiconductor junction: weak Fermi level pinning enables effective tuning of Schottky barrier}, Sci.\ Adv.\ \textbf{2}, e1600069 (2016), \doi{10.1126/sciadv.1600069}.",
"refliu18": r"Y.~Liu \textit{et al.}, \textit{Approaching the Schottky--Mott limit in van der Waals metal--semiconductor junctions}, Nature \textbf{557}, 696 (2018), \doi{10.1038/s41586-018-0129-8}.",
"refwangy19": r"Y.~Wang \textit{et al.}, \textit{Van der Waals contacts between three-dimensional metals and two-dimensional semiconductors}, Nature \textbf{568}, 70 (2019), \doi{10.1038/s41586-019-1052-3}.",
"refshen": r"P.-C.~Shen \textit{et al.}, \textit{Ultralow contact resistance between semimetal and monolayer semiconductors}, Nature \textbf{593}, 211 (2021), \doi{10.1038/s41586-021-03472-9}.",
"refli23": r"W.~Li \textit{et al.}, \textit{Approaching the quantum limit in two-dimensional semiconductor contacts}, Nature \textbf{613}, 274 (2023), \doi{10.1038/s41586-022-05431-4}.",
"refwangy22": r"Y.~Wang \textit{et al.}, \textit{P-type electrical contacts for 2D transition-metal dichalcogenides}, Nature \textbf{610}, 61 (2022), \doi{10.1038/s41586-022-05134-w}.",
"refchen": r"Y.-H.~Chen \textit{et al.}, \textit{Oxidized-monolayer tunneling barrier for strong Fermi-level depinning in layered InSe transistors}, npj 2D Mater.\ Appl.\ \textbf{3}, 49 (2019), \doi{10.1038/s41699-019-0133-3}.",
"refguo": r"Y.~Guo and J.~Robertson, \textit{Band structure, band offsets, substitutional doping, and Schottky barriers of bulk and monolayer InSe}, Phys.\ Rev.\ Mater.\ \textbf{1}, 044004 (2017), \doi{10.1103/PhysRevMaterials.1.044004}.",
"refshi": r"B.~Shi \textit{et al.}, \textit{n-Type ohmic contact and p-type Schottky contact of monolayer InSe transistors}, Phys.\ Chem.\ Chem.\ Phys.\ \textbf{20}, 24641 (2018), \doi{10.1039/C8CP04615H}.",
"refzhao": r"Q.~Zhao, W.~Jie, T.~Wang, A.~Castellanos-Gomez and R.~Frisenda, \textit{InSe Schottky diodes based on van der Waals contacts}, Adv.\ Funct.\ Mater.\ \textbf{30}, 2001307 (2020), \doi{10.1002/adfm.202001307}.",
"refwang": r"L.~Wang \textit{et al.}, Science \textbf{342}, 614 (2013), \doi{10.1126/science.1244358}.",
"refguim": r"M.H.D.~Guimar\~aes \textit{et al.}, ACS Nano \textbf{10}, 6392 (2016), \doi{10.1021/acsnano.6b02879}.",
"refcheng": r"Z.~Cheng, Y.~Yu, S.~Singh, K.~Price, S.G.~Noyce, Y.-C.~Lin, L.~Cao and A.D.~Franklin, \textit{Immunity to contact scaling in MoS$_2$ transistors using in situ edge contacts}, Nano Lett.\ \textbf{19}, 5077 (2019), \doi{10.1021/acs.nanolett.9b01355}.",
"refjain": r"A.~Jain \textit{et al.}, \textit{One-dimensional edge contacts to a monolayer semiconductor}, Nano Lett.\ \textbf{19}, 6914 (2019), \doi{10.1021/acs.nanolett.9b02166}.",
"refszabo": r"\'A.~Szab\'o, A.~Jain, M.~Parzefall, L.~Novotny and M.~Luisier, \textit{Electron transport through metal/MoS$_2$ interfaces: edge- or area-dependent process?}, Nano Lett.\ \textbf{19}, 3641 (2019), \doi{10.1021/acs.nanolett.9b00678}.",
"refparto": r"K.~Parto, A.~Pal, T.~Chavan, K.~Agashiwala, C.-H.~Yeh, W.~Cao and K.~Banerjee, \textit{One-dimensional edge contacts to two-dimensional transition-metal dichalcogenides: uncovering the role of Schottky-barrier anisotropy in charge transport across MoS$_2$/metal interfaces}, Phys.\ Rev.\ Appl.\ \textbf{15}, 064068 (2021), \doi{10.1103/PhysRevApplied.15.064068}.",
"refarutchelvan": r"G.~Arutchelvan \textit{et al.}, \textit{From the metal to the channel: a study of carrier injection through the metal/2D MoS$_2$ interface}, Nanoscale \textbf{9}, 10869 (2017), \doi{10.1039/C7NR02487H}.",
"refdasapp": r"S.~Das and J.~Appenzeller, \textit{Where does the current flow in two-dimensional layered systems?}, Nano Lett.\ \textbf{13}, 3396 (2013), \doi{10.1021/nl401831u}.",
"refdas13": r"S.~Das, H.-Y.~Chen, A.V.~Penumatcha and J.~Appenzeller, \textit{High performance multilayer MoS$_2$ transistors with scandium contacts}, Nano Lett.\ \textbf{13}, 100 (2013), \doi{10.1021/nl303583v}.",
"refpenu": r"A.V.~Penumatcha, R.B.~Salazar and J.~Appenzeller, \textit{Analysing black phosphorus transistors using an analytic Schottky barrier MOSFET model}, Nat.\ Commun.\ \textbf{6}, 8948 (2015), \doi{10.1038/ncomms9948}.",
"refsangwan": r"V.K.~Sangwan \textit{et al.}, \textit{Thickness-dependent charge transport in exfoliated indium selenide vertical field-effect transistors}, Appl.\ Phys.\ Lett.\ \textbf{115}, 243104 (2019), \doi{10.1063/1.5128808}.",
"refrai": r"A.~Rai, V.K.~Sangwan, J.T.~Gish, M.C.~Hersam and D.G.~Cahill, \textit{Anisotropic thermal conductivity of layered indium selenide}, Appl.\ Phys.\ Lett.\ \textbf{118}, 073101 (2021), \doi{10.1063/5.0042091}.",
"refyalon": r"E.~Yalon \textit{et al.}, \textit{Energy dissipation in monolayer MoS$_2$ electronics}, Nano Lett.\ \textbf{17}, 3429 (2017), \doi{10.1021/acs.nanolett.7b00252}.",
"reffiori": r"G.~Fiori \textit{et al.}, \textit{Electronics based on two-dimensional materials}, Nat.\ Nanotechnol.\ \textbf{9}, 768 (2014), \doi{10.1038/nnano.2014.207}.",
"refmudd": r"G.W.~Mudd \textit{et al.}, \textit{The direct-to-indirect band gap crossover in two-dimensional van der Waals indium selenide crystals}, Sci.\ Rep.\ \textbf{6}, 39619 (2016), \doi{10.1038/srep39619}.",
"refmudd2013": r"G.W.~Mudd \textit{et al.}, \textit{Tuning the bandgap of exfoliated InSe nanosheets by quantum confinement}, Adv.\ Mater.\ \textbf{25}, 5714 (2013), \doi{10.1002/adma.201302616}.",
"refmag": r"S.J.~Magorrian, V.~Z\'olyomi and V.I.~Fal'ko, \textit{Electronic and optical properties of two-dimensional InSe from a DFT-parametrized tight-binding model}, Phys.\ Rev.\ B \textbf{94}, 245431 (2016), \doi{10.1103/PhysRevB.94.245431}.",
"refbednarczyk": r"D.~Bednarczyk and J.~Bednarczyk, \textit{The approximation of the Fermi--Dirac integral $F_{1/2}(\eta)$}, Phys.\ Lett.\ A \textbf{64}, 409 (1978).",
"refcrowell": r"C.R.~Crowell and S.M.~Sze, \textit{Current transport in metal-semiconductor barriers}, Solid-State Electron.\ \textbf{9}, 1035 (1966), \doi{10.1016/0038-1101(66)90127-4}.",
"refieong": r"M.~Ieong, P.M.~Solomon, S.E.~Laux, H.-S.P.~Wong and D.~Chidambarrao, \textit{Comparison of raised and Schottky source/drain MOSFETs using a novel tunneling contact model}, IEDM Tech.\ Dig., 733 (1998), \doi{10.1109/IEDM.1998.746461}.",
"refpad": r"F.A.~Padovani and R.~Stratton, Solid-State Electron.\ \textbf{9}, 695 (1966), \doi{10.1016/0038-1101(66)90097-9}.",
"refwachutka": r"G.K.~Wachutka, \textit{Rigorous thermodynamic treatment of heat generation and conduction in semiconductor device modeling}, IEEE Trans.\ Comput.-Aided Des.\ \textbf{9}, 1141 (1990), \doi{10.1109/43.62751}.",
"refenglish": r"C.D.~English, G.~Shine, V.E.~Dorgan, K.C.~Saraswat and E.~Pop, \textit{Improved contacts to MoS$_2$ transistors by ultra-high vacuum metal deposition}, Nano Lett.\ \textbf{16}, 3824 (2016), \doi{10.1021/acs.nanolett.6b01309}.",
}
f = sys.argv[1]
s = open(f).read()
body = s.split(r"\begin{thebibliography}")[0]
order = []
for grp in re.findall(r"\\cite\{([^}]*)\}", body):
    for k in grp.split(","):
        k = k.strip()
        if k not in order:
            order.append(k)
missing = [k for k in order if k not in B]
unused = [k for k in B if k not in order]
print("cited:", len(order), "missing:", missing, "unused:", unused)
assert not missing
bib = "\n".join(rf"\bibitem{{{k}}} {B[k]}" for k in order)
s = s.replace("BIBLIOGRAPHY", bib)
open(f, "w").write(s)
