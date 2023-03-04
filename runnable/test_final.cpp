#include "pb/det/det_sql.h"
#include "pb/det/det_exp.h"
#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/dlv.h"
#include "pb/det/lsr.h"
#include "pb/det/sr.h"
#include "pb/det/kd_tree.h"

#include "pb/core/checker.h"
#include "pb/core/dual.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"

#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"

using namespace pb;

bool is_finished = false;
double max_RAM = 0;
const int kRAMSleep = 10; // In ms

void inspectRAM(){
  is_finished = false;
  max_RAM = 0;
  while (!is_finished){
    max_RAM = max(max_RAM, currentRAM());
    std::this_thread::sleep_for(std::chrono::milliseconds(kRAMSleep));
  }
}

/******************************************/
void Z1(){
  string exp_name = "Z1";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  exp.q = 1;
  int R = 10;
  int cnt = 0;
  for (auto g : exp.g2){
    exp.g = g;
    if (g < 0.01) continue;
    for (auto S : exp.S3){
      // if (cnt <= 4) continue;
      // if (S <= 100000) continue;
      exp.S = S;
      for (int i = 1; i <= R; i ++){
        cnt ++;
        if (S <= 1000000) continue;
        // if (cnt <= 2) continue;
        // if (i <= 1) continue;
        // if (S == 100000) continue;
        // if (S == 1000000 && i <= 3) continue;
        string label = fmt::format("LSR*_i{}_g{}_s{}", i, formatFloat(g), exp.S);
        exp.seed = i;
        cout << "det_sql generate\n";
        DetSql det_sql = exp.generate();
        cout << "dlv partition\n";
        double p_exe = exp.dlvPartition(false);
        exp.write(label + "_ptime", p_exe);
        cout << "lsr_prob generate\n";
        LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
        cout << "det_prob generate\n";
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        for (auto H : exp.H8){
          exp.H = H;
          cout << "lsr_bound_generate\n";
          lsr_prob.generateBounds(exp.Es[exp.q], exp.a, exp.H);
          cout << "det_prob_copy_bound\n";
          det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
          cout << "lsr\n";
          LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
          cout << label << " " << solMessage(lsr.status) << endl;
          if (lsr.status == Found){
            cout << "dual\n";
            Dual dual = Dual(exp.C, det_prob);
            // LsrChecker ch = LsrChecker(lsr_prob);
            // cout << feasMessage(ch.checkIlpFeasibility(lsr.ilp_sol)) << " " << feasMessage(ch.checkLpFeasibility(lsr.lp_sol)) << endl;
            // assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
            // assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
            exp.write(label + "_time", H, lsr.exe_ilp);
            exp.write(label + "_err", H, pctError(lsr.ilp_score, dual.score));
          }
        }
      }
    }
  }
}
/******************************************/
// void A1(){
//   string exp_name = "A1";
//   DetExp exp = DetExp(exp_name);
//   exp.o = 8;
//   int R = 1;
//   DetSql det_sql = exp.generate();
//   DetProb prob = DetProb(det_sql, -1, exp.seed);
//   for (auto C : exp.C7){
//     exp.C = C;
//     for (int i = 0; i < R; i ++){
//       prob.generateBounds(exp.E, exp.a, exp.H);
//       DualReducer dr = DualReducer(exp.C, prob);
//       exp.write("D", exp.C, dr.exe_lp);
//       exp.write("DR", exp.C, dr.exe_ilp);
//     }
//   }
//   for (auto C : exp.C7){
//     exp.C = C;
//     exp.partition_name = fmt::format("{}_C{}", exp_name, exp.C);
//     LsrProb lsr_prob = LsrProb(det_sql, exp.partition_name, exp.seed);
//     double p_exe = exp.partition(false);
//     exp.write("DLV", exp.C, p_exe);
//     for (int i = 0; i < R; i ++){
//       lsr_prob.generateBounds(exp.E, exp.a, exp.H);
//       LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
//       // cout << solMessage(lsr.status) << endl;
//       if (lsr.status == Found){
//         LsrChecker ch = LsrChecker(lsr_prob);
//         assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
//         assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
//         exp.write("LSR", exp.C, lsr.exe_ilp);
//       }
//     }
//   }
// }
// /******************************************/
// void A2(){
//   string exp_name = "A2";
//   DetExp exp = DetExp(exp_name);
//   int R = 3;
//   PgManager pg = PgManager();
//   for (auto o : exp.o6){
//     // if (o > 7) continue;
//     double N = pow(10, o);
//     exp.o = o;
//     exp.partition_name = fmt::format("{}_{}", exp_name, exp.o);
//     string kd_partition_name = "kd_" + exp.partition_name;
//     DetSql det_sql = exp.generate();
//     DetProb det_prob;
//     vector<LsrProb> lsr_probs;
//     double ground = 0;
//     // cout << "WTF1\n";
//     if (o <= 9){
//       LsrProb lsr_prob = LsrProb(det_sql, exp.partition_name, exp.seed);
//       for (int i = 0; i < R; i ++){
//         exp.H = 1;
//         lsr_prob.generateBounds(exp.E, exp.a, exp.H);
//         lsr_probs.push_back(lsr_prob);
//       }
//       if (o <= 7){
//         //sr partition
//         int tau = (int) (0.001 * N);
//         // cout << o << " " << tau << endl;
//         pg.dropTable("[1P]_" + exp.getTableName() + "_" + kd_partition_name); 
//         pg.dropTable("[1G]_" + exp.getTableName() + "_" + kd_partition_name); 
//         KDTree kt;
//         kt.partitionTable(exp.getTableName(), kd_partition_name, exp.getCols(), tau, DBL_MAX);
//       }
//       if (o >= 6 && o <= 9){
//         //lsr partition
//         exp.partition();
//       }
//     }
//     if (o <= 8){
//       det_prob = DetProb(det_sql, -1, exp.seed);
//       for (int i = 0; i < R; i ++){
//         if (o < 5) det_prob.generateBounds(exp.E, exp.a, exp.H);
//         else {
//           LsrProb lsr_prob = lsr_probs[i];
//           det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
//         }
//         if (o <= 7){
//           GurobiSolver gs = GurobiSolver(det_prob);
//           gs.solveLp();
//           exp.write("GD", N, gs.exe_lp);
//           ground = gs.lp_score;
//           exp.write("GD_aux", N, pctError(gs.lp_score, ground));
//         }
//         if (o <= 6 && 0==1){
//           GurobiSolver gs = GurobiSolver(det_prob);
//           gs.solveIlp(0);
//           exp.write("GDR", N, gs.exe_ilp);
//           exp.write("GDR_aux", N, pctError(gs.ilp_score, ground));
//         }
//         DualReducer dr = DualReducer(exp.C, det_prob);
//         if (dr.status == Found){
//           if (o == 8) ground = dr.lp_score;
//           exp.write("D", N, dr.exe_lp);
//           exp.write("D_aux", N, pctError(dr.lp_score, ground));
//           exp.write("DR", N, dr.exe_ilp);
//           exp.write("DR_aux", N, pctError(dr.ilp_score, ground));
//         }
//         if (o <= 7){
//           LsrProb lsr_prob = lsr_probs[i];
//           lsr_prob.partition_name = kd_partition_name;
//           SketchRefine sr = SketchRefine(lsr_prob);
//           map<long long, long long> sol;
//           bool is_success = sr.sketchAndRefine(sol);
//           LsrChecker ch = LsrChecker(lsr_prob);
//           if (is_success) exp.write("SR", N, sr.exec_sr);
//           if (is_success && ch.checkIlpFeasibility(sol)==Feasibility){
//             double score = ch.getScore(sol);
//             exp.write("SR_aux", N, pctError(score, ground));
//           }
//         }
//         if (o >= 6 && o <= 8){
//           LsrProb lsr_prob = lsr_probs[i];
//           lsr_prob.partition_name = exp.partition_name;
//           LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
//           if (lsr.status == Found){
//             LsrChecker ch = LsrChecker(lsr_prob);
//             assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
//             assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
//             exp.write("LSR", N, lsr.exe_ilp);
//             exp.write("LSR_aux", N, pctError(lsr.ilp_score, ground));
//           }
//         }
//       }
//     }
//     if (o == 9){
//       for (int i = 0; i < R; i ++){
//         LsrProb lsr_prob = lsr_probs[i];
//         LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
//         if (lsr.status == Found){
//           LsrChecker ch = LsrChecker(lsr_prob);
//           assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
//           assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
//           ground = lsr.lp_score;
//           exp.write("LSR", N, lsr.exe_ilp);
//           exp.write("LSR_aux", N, pctError(lsr.ilp_score, ground));
//         }
//       }
//     }
//   }
// }
// /******************************************/
void A3(){
  PgManager pg = PgManager();
  string exp_name = "A3";
  DetExp exp = DetExp(exp_name);
  exp.o = 6;
  int R = 10;
  // double N = pow(10, exp.o);
  for (int q = 0; q < 2; q ++){
    // if (q == 0) continue;
    exp.q = q;
    for (int i = 1; i <= R; i ++){
      // if (q == 0 && i == 21) continue;
      // if (q == 1 && i == 9) continue;
      // if (q == 0 && i <= 13) continue;
      // if (q == 0) continue;
      // if (q == 1 and i <= 9) continue;
      // cout << "Seed: " << i << endl;
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.kdPartition();
      exp.dlvPartition();
      for (auto H : exp.H8){
        // if (i == 14 && H <= 13) continue;
        exp.H = H;
        exp.E = exp.Es[exp.q];
        string label = fmt::format("{}_{}", exp.datasets[exp.q], exp.E);
        LsrProb lsr_prob = LsrProb(det_sql, "", exp.seed);
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
        lsr_prob.partition_name = exp.getDlvPartitionName();
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        cout << "LSR finished\n";
        lsr_prob.partition_name = exp.getKdPartitionName();
        SketchRefine sr = SketchRefine(lsr_prob);
        map<long long, long long> sr_sol;
        bool is_success = sr.sketchAndRefine(sr_sol);
        cout << "SR finished\n";
        LsrChecker ch = LsrChecker(lsr_prob);
        Checker det_ch = Checker(det_prob);
        double ground = 0;
        bool is_positive = false;
        // GurobiSolver gs_lp = GurobiSolver(det_prob);
        // gs_lp.solveLp();
        // cout << "GS LP finished\n";
        // ground = gs_lp.lp_score;
        // It is mathematically proved that dr.status or lsr.status == Found means Checker is correct
        
        // if (lsr.status == Found){
        //   is_positive = true;
        // } else{
        //   cout << "Start gs decide has solution or not\n"; 
        //   GurobiSolver gs = GurobiSolver(det_prob);
        //   is_positive = gs.hasIlpSolution();
        //   cout << "gs finished\n";
        // }
        cout << "Start dual\n";
        Dual dual = Dual(exp.C, det_prob);
        ground = dual.score;
        cout << "Start gs\n"; 
        GurobiSolver gs = GurobiSolver(det_prob);
        gs.solveIlp();
        if (gs.ilp_status == Found){
          is_positive = true;
          exp.write("GR_err_" + label, exp.H, pctError(gs.ilp_score, ground));
        } else {
          is_positive = false;
        }

        if (is_success && ch.checkIlpFeasibility(sr_sol)==Feasibility){
          double sr_score = ch.getScore(sr_sol);
          exp.write("SR_TP_" + label, exp.H, 1);
          exp.write("SR_err_" + label, exp.H, pctError(sr_score, ground));
        } else{
          if (is_positive) exp.write("SR_FN_" + label, exp.H, 1);
          else exp.write("SR_TN_" + label, exp.H, 1);
        }
        if (lsr.status == Found){
          exp.write("LSR_TP_" + label, exp.H, 1);
          exp.write("LSR_err_" + label, exp.H, pctError(lsr.ilp_score, ground));
        } else{
          if (is_positive) exp.write("LSR_FN_" + label, exp.H, 1);
          else exp.write("LSR_TN_" + label, exp.H, 1);
        }
      }
    }
  }
}
// /******************************************/
// void P1(){
//   string exp_name = "P1";
//   DetExp exp = DetExp(exp_name);
//   exp.o = 8;
//   DetSql det_sql = exp.generate();
//   for (double M : exp.M6){
//     exp.M = M;
//     std::thread ram (inspectRAM);
//     exp.partition_name = fmt::format("{}_M{}", exp_name, exp.M);
//     double p_exe = exp.partition(false);
//     is_finished = true;
//     ram.join();
//     exp.write("DLV", exp.M, max_RAM);
//     exp.write("DLV_aux", exp.M, p_exe);
//   }
// }
// /******************************************/
// void P2(){
//   string exp_name = "P2";
//   DetExp exp = DetExp(exp_name);
//   for (int o : exp.o4){
//     exp.o = o;
//     exp.partition_name = fmt::format("{}_o{}", exp_name, exp.o);
//     DynamicLowVariance dlv = DynamicLowVariance(exp.C, exp.g, exp.M, exp.S);
//     string table_name = exp.getTableName();
//     dlv.dropPartition(table_name, exp.partition_name);
//     dlv.partition(table_name, exp.partition_name);
//     double N = pow(10.0, o);
//     exp.write("DLV_read", N, dlv.pro.time(0));
//     exp.write("DLV_write", N, dlv.pro.time(1));
//     exp.write("DLV_index", N, dlv.pro.time(2));
//     exp.write("DLV_dlv", N, dlv.exe - dlv.pro.time(0) - dlv.pro.time(1) - dlv.pro.time(2));
//     exp.write("DLV_aux", N, dlv.exe);
//   }
// }
// /******************************************/
void L1(){
  string exp_name = "L1";
  DetExp exp = DetExp(exp_name);
  int R = 3;
  for (int o : exp.o4){
    exp.o = o;
    double N = pow(10.0, o);
    exp.H = 7;
    if (o <= 8) R = 5;
    else R = 3;
    for (int i = 0; i < R; i ++){
      exp.seed = i+1;
      exp.dlvPartition();
      DetSql det_sql = exp.generate();
      LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
      lsr_prob.generateBounds(exp.E, exp.a, exp.H);
      LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
      if (lsr.status == Found){
        lsr.pro.print();
        exp.write("LSR_read", N, lsr.pro.time(0));
        exp.write("LSR_dual", N, lsr.exe_dual);
        exp.write("LSR_gb", N, lsr.exe_gb);
        exp.write("LSR_sort", N, lsr.exe_sort);
        exp.write("LSR_misc", N, lsr.exe_ilp-lsr.exe_sort-lsr.exe_gb-lsr.exe_dual-lsr.pro.time(0));
        exp.write("LSR_total", N, lsr.exe_ilp);
      }
    }
  }
}
// /******************************************/
void L2(){
  string exp_name = "L2";
  DetExp exp = DetExp(exp_name);
  exp.o = 8;
  int R = 3;
  for (int q = 0; q < 2; q ++){
    if (q == 1) continue;
    exp.q = q;
    string label = "LSR_" + exp.datasets[q];
    for (int i = 0; i < R; i ++){
      exp.seed = i+1;
      exp.dlvPartition();
      for (auto F : exp.F5){
        // if (F != 0.7) continue;
        exp.F = F;
        DetSql det_sql = exp.generate();
        det_sql.addFilterWithRatio(exp.filtered_cols[exp.q], exp.F, Bounded);
        // cout << "OK-3\n";
        LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
        // cout << "OK-2\n";
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
        DualReducer dr = DualReducer(exp.C, det_prob);
        // cout << label << " " << solMessage(lsr.status) << endl;
        if (lsr.status == Found){
          LsrChecker ch = LsrChecker(lsr_prob);
          assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
          assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
          double ground = dr.lp_score;
          exp.write(label+"_lp", exp.F, pctError(lsr.lp_score, ground));
          exp.write(label+"_ilp", exp.F, pctError(lsr.ilp_score, ground));
          exp.write(label+"_aux", exp.F, lsr.exe_ilp);
        }
      }
    }
  }
}
// /******************************************/
void A4(){
  string exp_name = "A4";
  DetExp exp = DetExp(exp_name);
  PgManager pg = PgManager();
  for (int q = 0; q < 2; q ++){
    if (q == 0) continue;
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (auto H : exp.H4){
      exp.H = H;
      for (auto o : exp.o6){
        if (o == 9 && q == 1) continue;
        int R = 5;
        if (o == 9) R = 3;
        if (q == 0 && H <= 1) continue;
        if (q == 0 && H <= 3 && o <= 8) continue;
        // if (q == 0 && H == 7 && o <= 6) continue;
        exp.o = o;
        int seed_count = 0;
        exp.seed = 1;
        double N = pow(10, exp.o);
        string label = fmt::format("{}_H{}", exp.datasets[exp.q], exp.H);
        while (seed_count < R){
          DetSql det_sql = exp.generate();

          if (o <= 8){
            cout << "Start kd " + label << endl;
            exp.kdPartition(true);
          }

          ////////
          cout << "Start dlv "  + label << endl;
          exp.dlvPartition(true);
          ////////

          cout << "lsr_prob\n";
          LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
          lsr_prob.generateBounds(exp.E, exp.a, exp.H);
          DetProb det_prob;

          if (o <= 8){
            cout << "det_prob\n";
            det_prob = DetProb(det_sql, -1, exp.seed);
            det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
          }

          cout << "Start lsr " + label << endl;
          lsr_prob.partition_name = exp.getDlvPartitionName();
          LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
          if (lsr.status == Found){
            seed_count ++;
            LsrChecker ch = LsrChecker(lsr_prob);
            assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
            assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);

            double ground = 0;
            if (o <= 8){
              cout << "dual\n";
              Dual dual = Dual(exp.C, det_prob);
              ground = dual.score;
            } else ground = lsr.lp_score;

            if (o <= 6){
              cout << "Start gs " + label << endl;
              GurobiSolver gs = GurobiSolver(det_prob);
              gs.solveIlp(1e-4, kTimeLimit);
              if (gs.ilp_status == Found){
                exp.write("GDR_time_" + label, N, gs.exe_ilp);
                exp.write("GDR_pcterror_" + label, N, pctError(gs.ilp_score, ground));
              }
            }

            exp.write("LSR_time_" + label, N, lsr.exe_ilp);
            exp.write("LSR_pcterror_" + label, N, pctError(lsr.ilp_score, ground));

            if (o <= 8){
              cout << "Start sr " + label << endl;
              lsr_prob.partition_name = exp.getKdPartitionName();
              SketchRefine sr = SketchRefine(lsr_prob);
              map<long long, long long> sr_sol;
              bool is_success = sr.sketchAndRefine(sr_sol);
              LsrChecker ch = LsrChecker(lsr_prob);
              cout << "Sr finished " + label << endl; 
              if (is_success && ch.checkIlpFeasibility(sr_sol)==Feasibility){
                double sr_score = ch.getScore(sr_sol);
                exp.write("SR_time_" + label, N, sr.exec_sr);
                exp.write("SR_pcterror_" + label, N, pctError(sr_score, ground));
              }
            }

          }

          exp.seed ++;
        }
      }
    }
  }
}
/******************************************/
void test(){
  // DynamicLowVariance dlv = DynamicLowVariance();
  // auto ps = dlv.getPartitions();
  // for (auto pi : ps){
  //   string table_name = pi.first;
  //   auto partition_names = pi.second;
  //   for (auto partition_name : partition_names){
  //     dlv.dropPartition(table_name, partition_name);
  //   }
  // }

  // PgManager pg = PgManager();
  // auto tables = pg.listTables();
  // for (string table : tables){
  //   if (getLayerIndex(table) > 0){
  //     pg.dropTable(table);
  //   }
  // }

  PgManager pg = PgManager();
  auto tables = pg.listTables();
  for (string table : tables){
    if (pg.getSize(table) < 0.13e9) pg.dropTable(table);
  }

  // Inverse hardness test
  // string exp_name = "A4";
  // DetExp exp = DetExp(exp_name);
  // exp.o = 7;
  // DetSql det_sql = exp.generate();
  // DetProb det_prob = DetProb(det_sql, -1, exp.seed);
  // det_prob.generateBounds(exp.E, 0, 11.2);
  // double minE = 0;
  // double h = det_prob.det_bound.minHardness(minE, det_prob.bl, det_prob.bu);
  // cout << h << endl;
}
/******************************************/

int main() {
  // Z1();
  // A1();
  // A2();
  // A3();
  // A4();
  // P1();
  // P2();
  L1();
  // L2();
  // test();
  return 0;
}