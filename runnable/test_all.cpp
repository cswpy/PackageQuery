#include "pb/det/det_sql.h"
#include "pb/det/det_exp.h"
#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/dlv.h"
#include "pb/det/lsr.h"

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

void A1(){
  DetExp exp = DetExp("A1B");
  exp.o = 8;
  DetSql det_sql = exp.generate();
  DetProb prob = DetProb(det_sql, -1, exp.seed);
  for (auto C : exp.C6){
    for (auto H : exp.H3){
      exp.C = C;
      exp.H = H;
      prob.generateBounds(exp.E, exp.a, exp.H);
      DualReducer DR = DualReducer(exp.C, prob);
      exp.write("D", exp.C, DR.exe_lp);
      exp.write("DR", exp.C, DR.exe_ilp);
    }
  }
}

void A2(){
  DetExp exp = DetExp("A2");
  exp.o = 8;
  DetSql det_sql = exp.generate();
  for (auto N : exp.N5){
    DetProb prob = DetProb(det_sql, N, exp.seed);
    for (auto H : exp.H3){
      exp.H = H;
      prob.generateBounds(exp.E, exp.a, exp.H);
      DualReducer DR = DualReducer(exp.C, prob);
      exp.write("D", N, DR.exe_lp);
      exp.write("DR", N, DR.exe_ilp);
      if (DR.status == Found){
        double ground = DR.lp_score;
        exp.write("D_aux", N, pctError(DR.lp_score, ground));
        exp.write("DR_aux", N, pctError(DR.ilp_score, ground));
        if (N <= 10000000){
          GurobiSolver gs = GurobiSolver(prob);
          gs.solveLp();
          exp.write("GLP", N, gs.exe_lp);
          exp.write("GLP_aux", N, pctError(gs.lp_score, ground));
        }
        if (N <= 1000000){
          GurobiSolver gs = GurobiSolver(prob);
          gs.solveIlp();
          exp.write("GILP", N, gs.exe_ilp);
          exp.write("GILP_aux", N, pctError(gs.ilp_score, ground));
        }
      }
    }
  }
}

void A3(){
  DetExp exp = DetExp("A3");
  exp.o = 7;
  long long N = 100000;
  int R = 10;
  for (int seed = 1; seed <= R; seed ++){
    for (int q = 0; q < 2; q ++){
      exp.q = q;
      DetSql det_sql = exp.generate();
      DetProb prob = DetProb(det_sql, N, seed);
      for (auto E : exp.E2){
        for (auto H : exp.H8){
          exp.H = H;
          exp.E = E;
          prob.generateBounds(exp.E, exp.a, exp.H);
          DualReducer DR = DualReducer(exp.C, prob);
          double ground = DR.lp_score;
          int relevant = 1;
          if (DR.status == Found) exp.write(fmt::format("DR_{}_{}", exp.q, exp.E), exp.H, pctError(DR.ilp_score, ground));
          else{
            GurobiSolver gs = GurobiSolver(prob);
            if (gs.hasIlpSolution(10.0)) relevant = 0;
          }
          exp.write(fmt::format("DR_{}_{}_aux", exp.q, exp.E), exp.H, relevant);
        }
      }
    }
  }
}

/******************************************/

void P1(){
  DetExp exp = DetExp("P1");
  exp.o = 8;
  for (double M : exp.M6){
    exp.M = M;
    std::thread ram (inspectRAM);
    DynamicLowVariance dlv = DynamicLowVariance(exp.C, exp.g, exp.M);
    dlv.partition(exp.getTableName(), "M" + to_string(M));
    is_finished = true;
    ram.join();
    exp.write("DLV", exp.M, dlv.exe);
    exp.write("DLV_aux", exp.M, max_RAM);
  }
}

void P2(){
  DetExp exp = DetExp("P2");
  exp.o = 7;
  for (int C : exp.C6){
    exp.C = C;
    DynamicLowVariance dlv = DynamicLowVariance(exp.C, exp.g, exp.M);
    dlv.partition(exp.getTableName(), "C" + to_string(exp.C));
    exp.write("DLV", exp.C, dlv.exe);
  }
}

void P3(){
  DetExp exp = DetExp("P3");
  for (int o : exp.o4){
    exp.o = o;
    DynamicLowVariance dlv = DynamicLowVariance(exp.C, exp.g, exp.M);
    dlv.partition(exp.getTableName(), "P3");
    double X = pow(10.0, o);
    exp.write("DLV_read", X, dlv.pro.time(0));
    exp.write("DLV_write", X, dlv.pro.time(1));
    exp.write("DLV_index", X, dlv.pro.time(2));
    exp.write("DLV_aux", X, dlv.exe);
  }
}

/******************************************/

void L1(){
  DetExp exp = DetExp("L1");
  for (int o : exp.o4){
    exp.o = o;
    exp.partition_name = "P3";
    DetSql det_sql = exp.generate();
    double X = pow(10.0, o);
    LsrProb lsr_prob = LsrProb(det_sql, exp.partition_name, exp.seed);
    for (auto H : exp.H3){
      exp.H = H;
      lsr_prob.generateBounds(exp.E, exp.a, exp.H);
      LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
      if (lsr.status == Found){
        LsrChecker ch = LsrChecker(lsr_prob);
        assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
        assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
        exp.write("LSR", X, lsr.exe_ilp);
        double ground = lsr.lp_score;
        if (o <= 8){
          DetProb det_prob = DetProb(det_sql, -1, exp.seed);
          det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
          DualReducer dr = DualReducer(exp.C, det_prob);
          ground = dr.lp_score;
          exp.write("DR", X, dr.exe_ilp);
          exp.write("DR_aux", X, pctError(dr.ilp_score, ground));
        }
        exp.write("LSR_aux", X, pctError(lsr.ilp_score, ground));
      }
    }
  }
}

void L2(){
  DetExp exp = DetExp("L2");
  exp.o = 8;
  for (int C : exp.C6){
    exp.C = C;
    exp.partition_name = "P3";
    DetSql det_sql = exp.generate();
    LsrProb lsr_prob = LsrProb(det_sql, exp.partition_name, exp.seed);
    for (auto H : exp.H3){
      exp.H = H;
      lsr_prob.generateBounds(exp.E, exp.a, exp.H);
      LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
      if (lsr.status == Found){
        LsrChecker ch = LsrChecker(lsr_prob);
        assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
        assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
        exp.write("LSR", exp.C, lsr.exe_ilp);
      }
    }
  }
}

void L3(){
  DetExp exp = DetExp("L3");
  exp.o = 7;
  int R = 10;
  for (int seed = 1; seed <= R; seed ++){
    for (int q = 0; q < 2; q ++){
      exp.q = q;
      exp.seed = seed;
      DetSql det_sql = exp.generate();
      exp.partition_name = "L3x" + to_string(seed);
      exp.partition();
      for (auto E : exp.E2){
        for (double H : exp.H8){
          exp.E = E;
          exp.H = H;
          LsrProb lsr_prob = LsrProb(det_sql, exp.partition_name, exp.seed);
          lsr_prob.generateBounds(exp.E, exp.a, exp.H);
          LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
          double ground = lsr.lp_score;
          int relevant = 1;
          if (lsr.status == Found){
            LsrChecker ch = LsrChecker(lsr_prob);
            assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
            assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
            exp.write(fmt::format("LSR_{}_{}", exp.q, exp.E), exp.H, pctError(lsr.ilp_score, ground));
          } else{
            DetProb det_prob = DetProb(det_sql, -1, exp.seed);
            det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
            GurobiSolver gs = GurobiSolver(det_prob);
            if (gs.hasIlpSolution(10.0)) relevant = 0;
          }
          exp.write(fmt::format("LSR_{}_{}_aux", exp.q, exp.E), exp.H, relevant);
        }
      }
    }
  }
}

void L4(){
  DetExp exp = DetExp("L4");
  for (int o : exp.o4){
    exp.o = o;
    exp.partition_name = "P3";
    DetSql det_sql = exp.generate();
    double X = pow(10.0, o);
    LsrProb lsr_prob = LsrProb(det_sql, exp.partition_name, exp.seed);
    for (auto H : exp.H3){
      exp.H = H;
      lsr_prob.generateBounds(exp.E, exp.a, exp.H);
      LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
      if (lsr.status == Found){
        LsrChecker ch = LsrChecker(lsr_prob);
        assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
        assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
        exp.write("LSR_read", X, lsr.pro.time(0));
        exp.write("LSR_aux", X, lsr.exe_ilp);
      }
    }
  }
}

/******************************************/

int main(){
  A1();
  // A2();
  // A3();
  // P1();
  // P2();
  // P3();
  // L1();
  // L2();
  // L3();
  L4();
}