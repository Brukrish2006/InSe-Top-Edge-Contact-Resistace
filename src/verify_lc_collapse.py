import devsim
import csv
from sclc_test import run_simulation

def main():
    wf_52 = 5.2
    mu = 50.0
    kappa = 0.085
    thickness = 5e-7
    L_gap_50 = 50e-7
    L_c_50 = 50e-7
    L_c_10 = 10e-7
    
    print("Running Top Lc=50nm mesh=2.0...")
    devsim.reset_devsim()
    v, _, r_top_50 = run_simulation('dev_top50_m2', 'top', L_c_50, L_gap_50, thickness, wf_52, mu, kappa, 2.0, 0.1, 0.1, 0.1)
    
    print("Running Top Lc=10nm mesh=2.0...")
    devsim.reset_devsim()
    _, _, r_top_10 = run_simulation('dev_top10_m2', 'top', L_c_10, L_gap_50, thickness, wf_52, mu, kappa, 2.0, 0.1, 0.1, 0.1)
    
    print("Running Edge mesh=2.0...")
    devsim.reset_devsim()
    _, _, r_edge = run_simulation('dev_edge_m2', 'edge', thickness, L_gap_50, thickness, wf_52, mu, kappa, 2.0, 0.1, 0.1, 0.1)
    
    with open('verify_lc_collapse.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['V_DS', 'R_Top_Lc50_m2', 'R_Top_Lc10_m2', 'R_Edge_m2'])
        writer.writerow([v[0], r_top_50[0], r_top_10[0], r_edge[0]])
        
    print(f"Results at 0.1V: R_Top50={r_top_50[0]}, R_Top10={r_top_10[0]}, R_Edge={r_edge[0]}")

if __name__ == '__main__':
    main()
