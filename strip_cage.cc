#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TLine.h>
#include <TColor.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include "TSystem.h"
#include "Garfield/ComponentElmer.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/Random.hh"
#include "Garfield/MediumGas.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/DriftLineRKF.hh"
#include <csignal>
#include <csetjmp>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TStyle.h>
#include <TPad.h>

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace Garfield;



void PlotFieldPointMap(Garfield::Sensor& sensor,
                       const std::string& out_dir,
                       const std::string& tag,
                       double d    = 0.1,      // characteristic spacing [cm]
                       double xlo  = -2.0, double xhi =  2.0,
                       double ylo  = -2.0, double yhi =  2.0,
                       double zmin = 0.0,  double zmax = 0.0,   // 0,0 = autoscale
                       double zprobe = 0.005,
                       int canvasW = 700) {

  const double r = 0.375 * d;                 // radius -> diameter = (3/4)*d

  std::vector<double> px, py, pe;
  double emin =  1e30, emax = -1e30;




  for (double x = xlo + 0.5 * d; x <= xhi; x += d) {
    for (double y = ylo + 0.5 * d; y <= yhi; y += d) {
      double ex = 0., ey = 0., ez = 0., v = 0.;
      Garfield::Medium* med = nullptr;
      int status = 0;
      sensor.ElectricField(x, y, zprobe, ex, ey, ez, v, med, status);
      if (status != 0 || med == nullptr) continue;   // inside metal or off-mesh
      const double e = std::sqrt(ex * ex + ey * ey + ez * ez);
      px.push_back(x); py.push_back(y); pe.push_back(e);
      if (e < emin) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (px.empty()) {
    std::cout << "PlotFieldPointMap: no valid points sampled.\n";
    return;
  }
  if (zmax > zmin) { emin = zmin; emax = zmax; }        // manual clip
  if (emax <= emin) emax = emin + 1.0;


  std::cout << "  field point map: " << px.size() << " points, |E| range "
            << emin << " .. " << emax << " V/cm\n";



  gStyle->SetPalette(kBird);
  const int nCol = TColor::GetNumberOfColors();

  // keep the dots round
  // canvas aspect follows the data aspect
  int canvasH = int(canvasW * (yhi - ylo) / (xhi - xlo));
  if (canvasH < 300)  canvasH = 300;
  if (canvasH > 2400) canvasH = 2400;

  TCanvas* c = new TCanvas(("cFieldPts_" + tag).c_str(),
                           "E field point map", canvasW, canvasH);
  c->SetRightMargin(0.18);
  c->SetLeftMargin(0.14);

  // supplies the axes and the palette axis on the right.
  TH2D* hF = new TH2D(("hFieldFrame_" + tag).c_str(),
                      ";x [cm];y [cm]", 10, xlo, xhi, 10, ylo, yhi);
  hF->SetStats(0);
  hF->SetBinContent(1, 1, emin);      // one entry so ROOT draws the palette
  hF->SetMinimum(emin);
  hF->SetMaximum(emax);
  hF->GetZaxis()->SetTitle("|E| [V/cm]");
  hF->Draw("COLZ");
  gPad->Update();

  // wipe the dummy bin, keep frame + palette
  TBox* blank = new TBox(xlo, ylo, xhi, yhi);
  blank->SetFillColor(kWhite);
  blank->Draw();

  for (size_t i = 0; i < px.size(); ++i) {
    int ic = int(nCol * (pe[i] - emin) / (emax - emin));
    if (ic < 0)      ic = 0;
    if (ic >= nCol)  ic = nCol - 1;
    const int col = TColor::GetColorPalette(ic);
    TEllipse* el = new TEllipse(px[i], py[i], r, r);
    el->SetFillColor(col);
    el->SetLineColor(col);
    el->Draw();
  }


  // strip planes for reference
  TLine* sL = new TLine(-1.5, ylo, -1.5, yhi);
  TLine* sR = new TLine( 1.5, ylo,  1.5, yhi);
  sL->SetLineColor(kGray + 2); sL->SetLineWidth(2); sL->Draw();
  sR->SetLineColor(kGray + 2); sR->SetLineWidth(2); sR->Draw();


  gPad->RedrawAxis();
  c->SaveAs((out_dir + "/strip_field_point_map_" + tag + ".pdf").c_str());
  delete c;
}






int main(int argc, char * argv[]) {
    TApplication app("app", &argc, argv);
    gROOT->SetBatch(true);








    // Create timestamped output directory
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time), "%Y-%m-%d_%H-%M-%S");
    std::string out_dir = "Outputs/" + ss.str();
   
    gSystem->MakeDirectory("Outputs");
    gSystem->MakeDirectory(out_dir.c_str());
    std::cout << "\nAll outputs will be saved to: " << out_dir << "\n";
#ifdef _OPENMP
    std::cout << "OpenMP enabled: up to " << omp_get_max_threads()
              << " threads (override with OMP_NUM_THREADS)\n";
#else
    std::cout << "OpenMP NOT enabled - running single-threaded\n";
#endif

    MediumGas gas;

    ComponentElmer elm("strip_cage/mesh.header", "strip_cage/mesh.elements", "strip_cage/mesh.nodes",
                       "dielectrics.dat", "strip_cage/strip_cage.result", "cm");
    elm.SetMedium(0, &gas);
    elm.EnablePeriodicityZ();


    Sensor sensor;
    sensor.AddComponent(&elm);
   
    // The strips now have rounded corners, eliminating the mathematical
    // singularities. We can now safely set the sensor area all the way out to x = +/- 2.0,
    // exactly matching field_cage.cc, allowing electrons to drift safely behind the strips!
    sensor.SetArea(-2.0, -16.0, -100.0, 2.0, 16.0, 100);

    AvalancheMC drift;
    drift.SetSensor(&sensor);




    drift.SetDistanceSteps(0.02);




    drift.SetTimeWindow(0.0, 40000.0);
   
    ViewDrift driftView;
    driftView.SetArea(-2.0, -16.0, -1.0, 2.0, 16.0, 1.0);



    // ======================================================================
    // SWITCH: 0 = Visual Diagnostic (67/33), 1 = Multi-Gas Spatial Sweep
    // ======================================================================
    int switch_mode = 1;

    if (switch_mode == 0) {
        std::cout << "\n=========================================\n";
        std::cout << "MODE 0: VISUAL DIAGNOSTIC (67% CF4)\n";
        std::cout << "=========================================\n";




        const std::string gas_file_m0 = "gas_tables/cf4_80_iso_20.gas";
        if (!gas.LoadGasFile(gas_file_m0)) {
            std::cerr << "FATAL: could not load " << gas_file_m0 << "\n";
            gSystem->Exit(1);
        }
        std::cout << "-> Gas: " << gas_file_m0 << "\n";
        drift.EnablePlotting(&driftView);




        const int nElectrons = 100;
        int leaked = 0;
       
        std::cout << "Drifting " << nElectrons << " random electrons for visualization...\n";
        for (int e = 0; e < nElectrons; ++e) {
            double x0 = -1.4 + RndmUniform() * 2.8;
            double y0 = -15 + RndmUniform() * 30;
            double z0 = 0.005, t0 = 0.0;


            std::cout << "  [e=" << e << "] x0=" << x0 << " y0=" << y0 << std::flush;


   
                        drift.DriftElectron(x0, y0, z0, t0);
                        const int nPoints = drift.GetNumberOfDriftLinePoints();
                        if (nPoints > 0) {
                            double x1, y1, z1, t1;
                            drift.GetDriftLinePoint(nPoints - 1, x1, y1, z1, t1);
                            if ((y1 < -15.99) && (std::abs(x1) < 1.5)) leaked++;
                            std::cout << " end (" << x1 << ", " << y1 << ") t=" << t1 << " ns, npts=" << nPoints << "\n";


                        }
                        std::cout << " -> ok\n";
                   



           
        }
        std::cout << "RESULTS: " << leaked << " out of " << nElectrons << " leaked past the field cage!\n\n";



        // 1. Drift Lines
        TCanvas* c1 = new TCanvas("c1", "Electron Drift Lines", 600, 600);
        driftView.SetCanvas(c1);
        driftView.Plot(true);
        TLine* stripL = new TLine(-1.5, -16.0, -1.5, 16.0);
        TLine* stripR = new TLine( 1.5, -16.0,  1.5, 16.0);
        stripL->SetLineColor(kGray+2); stripL->SetLineWidth(2); stripL->Draw();
        stripR->SetLineColor(kGray+2); stripR->SetLineWidth(2); stripR->Draw();
        c1->SaveAs((out_dir + "/strip_drift_lines_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/strip_drift_lines_67_33.pdf\n";




        // 2. Field Map
        TCanvas* cMap = new TCanvas("cMap", "Electric Field Contours", 800, 800);
        ViewField fieldView;
        fieldView.SetSensor(&sensor);
        fieldView.SetCanvas(cMap);
        fieldView.SetPlane(0, 0, 1, 0, 0, 0.005);
        fieldView.SetArea(-2.0, -16.0, 2.0, 16.0);
        fieldView.SetNumberOfContours(50);
        fieldView.SetNumberOfSamples2d(400, 400);
        fieldView.PlotContour("v");
        TLine* stripL2 = new TLine(-1.5, -16.0, -1.5, 16.0);
        TLine* stripR2 = new TLine( 1.5, -16.0,  1.5, 16.0);
        stripL2->SetLineColor(kGray+2); stripL2->SetLineWidth(2); stripL2->Draw();
        stripR2->SetLineColor(kGray+2); stripR2->SetLineWidth(2); stripR2->Draw();
        cMap->SaveAs((out_dir + "/strip_field_map_contours_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/strip_field_map_contours_67_33.pdf\n";




            // Zoomed: drift region near the middle, clipped colour range so the
        // bulk uniformity is actually visible.
        PlotFieldPointMap(sensor, out_dir, "zoom", 0.05,
                        -1.9, 1.9, -2.0, 2.0, 25.0, 45.0);




        // Top of the cage -- this is the region driving the efficiency upturn.
        PlotFieldPointMap(sensor, out_dir, "topend", 0.05,
                        -1.9, 1.9, 12.0, 16.0, 25.0, 45.0);




        // Full height, coarse grid, tall canvas.
        PlotFieldPointMap(sensor, out_dir, "full", 0.25,
                        -1.9, 1.9, -15.9, 15.9, 25.0, 45.0, 0.005, 320);





        // 3. Drift Velocity Graph
        ViewMedium mediumView;
        mediumView.SetMedium(&gas);
        TCanvas* cVel = new TCanvas("cVel", "Drift Velocity", 800, 600);
        mediumView.SetCanvas(cVel);
        mediumView.PlotElectronVelocity('e');
        cVel->SaveAs((out_dir + "/strip_drift_velocity_67_33.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/strip_drift_velocity_67_33.pdf\n";








        // 4. Mirror Symmetry Diagnostic
        // 4. Mirror Symmetry Diagnostic — scan the full drift volume
        double sum_diff = 0.0, sum_mag = 0.0;
        int nPts = 0;
        for (double yy = -15.9; yy <= 15.9; yy += 0.2) {
            for (double xx = 0.05; xx <= 1.49; xx += 0.05) {
                double ex1, ey1, ez1, v1, ex2, ey2, ez2, v2;
                Medium *m1 = nullptr, *m2 = nullptr;
                int s1 = 0, s2 = 0;
                elm.ElectricField( xx, yy, 0.005, ex1, ey1, ez1, v1, m1, s1);
                elm.ElectricField(-xx, yy, 0.005, ex2, ey2, ez2, v2, m2, s2);
                if (s1 != 0 || s2 != 0) continue;      // skip points outside the mesh
                sum_diff += std::abs(v1 - v2);
                sum_mag  += (std::abs(v1) + std::abs(v2)) / 2.0;
                ++nPts;
            }
        }
        std::cout << "-> Mirror Symmetry Diagnostic: "
                  << 100.0 * (1.0 - (sum_diff / sum_mag))
                  << "%  (" << nPts << " point pairs)\n";








        // ==========================================================
        // 5. Electric Field Uniformity Diagnostic
        // ==========================================================
        std::cout << "\n--- FIELD UNIFORMITY DIAGNOSTIC ---\n";
       
        double sum_e = 0.0;
        double sum_e2 = 0.0;
        int n_uniform_pts = 0;
       
        // Scan parameters
        double x_step = 0.1; // Check a column every 1 mm
        double y_step = .5; // Check down the Y axis every 2 mm
        double z_val = 0.005;



        for (double x = -1.3; x <= 1.301; x += x_step) {
            for (double y = -15.0; y <= 15.0; y += y_step) {
               
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
               
                elm.ElectricField(x, y, z_val, ex, ey, ez, v, m, status);
               
                if (status == 0) { // status == 0 means the point is safely inside the gas mesh
                    double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
                    sum_e += emag;
                    sum_e2 += (emag * emag);
                    n_uniform_pts++;
                }
            }
        }


        if (n_uniform_pts > 1) {
            double avg_e = sum_e / n_uniform_pts;
           
            // Standard deviation formula: sqrt( (Sum(x^2) - (Sum(x)^2 / N)) / (N - 1) )
            double variance = (sum_e2 - ((sum_e * sum_e) / n_uniform_pts)) / (n_uniform_pts - 1);
            if (variance < 0.0) variance = 0.0; // Prevent NaN from floating point rounding errors
            double std_dev = std::sqrt(variance);








            std::cout << "Sampled " << n_uniform_pts << " valid points inside X:[-1.49, 1.49], Y:[-15.0, 15.0]\n";
            std::cout << "Average Field (|E|) : " << avg_e << " V/cm\n";
            std::cout << "Standard Deviation  : " << std_dev << " V/cm\n";
           
            // Optional: Print it as a percentage of the average field (Coefficient of Variation)
            std::cout << "Field Fluctuation   : " << (std_dev / avg_e) * 100.0 << " %\n";
        } else {
            std::cout << "Error: Not enough valid mesh points found in the scanned region.\n";
        }
        std::cout << "-----------------------------------\n\n";








     
        {
        double ex = 0., ey = 0., ez = 0., v = 0.;
        Medium* m = nullptr;
        int status = 0;
        elm.ElectricField(1.0, -15.0, 0.005, ex, ey, ez, v, m, status);
        std::cout << "EDGE PROBE: |E|=" << sqrt(ex*ex+ey*ey+ez*ez) << " V/cm status=" << status << "\n";
        }
















    } else if (switch_mode == 1) {
        std::cout << "\n=========================================\n";
        std::cout << "MODE 1:  MULTI-GAS SWEEP\n";
        std::cout << "=========================================\n";


        // ==========================================================
        // 1. CALCULATE AVERAGE ELECTRIC FIELD
        // ==========================================================
        double sum_e = 0.0;
        int n_uniform_pts = 0;
        for (double x = -1.490; x <= 1.491; x += 0.1) {
            for (double y = -15.0; y <= 15.0; y += 0.2) {
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
                elm.ElectricField(x, y, 0.005, ex, ey, ez, v, m, status);
               
                if (status == 0) { // status == 0 means valid gas region
                    sum_e += std::sqrt(ex*ex + ey*ey + ez*ez);
                    n_uniform_pts++;
                }
            }
        }
        // Fallback to 33.0 if no points found, otherwise use calculated average
        double avg_e = (n_uniform_pts > 0) ? (sum_e / n_uniform_pts) : 33.0;
        std::cout << "-> Calculated Average Field for Sweep: " << avg_e << " V/cm\n";
















        // ==========================================================
        // 2. CONFIGURE GRAPHS
        // ==========================================================
        std::vector<int> cf4_fractions = {0, 20};
        double y_step = .5;
        // Sweep extent. Change these and the loop below follows.
        const double Y_MIN = -15.0;
        const double Y_MAX =  15.0;
        int colors[] = {kBlack, kBlue, kGreen+2, kOrange+1, kRed};








        // MultiGraphs for Spatial Efficiency
        TMultiGraph* mg = new TMultiGraph();
        mg->SetTitle("Total Spatial Efficiency vs Spawn Y;Spawn Y Position [cm];Electrons Passed [%]");
       
        TMultiGraph* mgCurves = new TMultiGraph();
        mgCurves->SetTitle("Spatial Efficiency Fits;Spawn Y Position [cm];Electrons Passed [%]");
       
        auto legend = new TLegend(0.80, 0.52, 0.985, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.035);








        // Transverse Diffusion Graph (Using dynamic average field)
        TGraph* grDiff = new TGraph();
        grDiff->SetTitle(Form("Transverse Diffusion at %.1f V/cm;CF4 Concentration [%%];Transverse Diffusion [cm^{1/2}]", avg_e));
        grDiff->SetMarkerStyle(20); grDiff->SetMarkerColor(kRed); grDiff->SetLineColor(kRed); grDiff->SetLineWidth(2);



        double eta = 0.0;
        gas.ElectronAttachment(0.0, avg_e, 0.0, 0., 0., 0., eta);
        std::cout << "attachment eta = " << eta << " /cm";
        if (eta > 1e-9) std::cout << "   (lambda = " << 1.0 / eta << " cm)";
        std::cout << "\n";






        // NEW: Drift Velocity Graph (Using dynamic average field)
        TGraph* grVel = new TGraph();
        grVel->SetTitle(Form("Drift velocity vs cf4/isobutane concentration at %.1f V/cm;CF4 Concentration [%%];Drift Velocity [cm/ns]", avg_e));
        grVel->SetMarkerStyle(20); grVel->SetMarkerColor(kMagenta+2); grVel->SetLineColor(kMagenta+2); grVel->SetLineWidth(2);








        int graph_idx = 0;








#ifndef _OPENMP
        // Serial path only. On the OpenMP path each thread builds its own
        // tracker inside the parallel region (see below).
        AvalancheMC local_drift;
        local_drift.SetSensor(&sensor);
        local_drift.SetDistanceSteps(0.02);
        local_drift.SetTimeWindow(0.0, 40000.0);   // ns — adjust to your real drift time
        local_drift.EnableAttachment();
#endif








        // ==========================================================
        // 3. RUN GAS SWEEP
        // ==========================================================
        for (int i : cf4_fractions) {
            std::string gas_file = "gas_tables/cf4_" + std::to_string(i) + "_iso_" + std::to_string(100 - i) + ".gas";
            // Checked, because a silent failure leaves the PREVIOUS mixture
            // loaded and produces a duplicate curve rather than an error.
            if (!gas.LoadGasFile(gas_file)) {
                std::cerr << "FATAL: could not load " << gas_file << "\n";
                gSystem->Exit(1);
            }








            std::string out_name = out_dir + "/epass_" + std::to_string(i) + "_" + std::to_string(100 - i) + ".txt";
            std::ofstream outFile(out_name);
            outFile << "Spawn_Y\tElectrons_Passed_Pct\n";
            std::string lost_name = out_dir + "/lost_" + std::to_string(i) + "_" + std::to_string(100 - i) + ".txt";
            std::ofstream lostFile(lost_name);
            lostFile << "Spawn_Y\tx0\tx_end\ty_end\tt_end\tnpts\n";
            std::cout << "\nTesting Gas: CF4 " << i << "% | Isobutane " << (100 - i) << "%\n";








            // Record Transverse Diffusion at avg_e
            double dl = 0, dt = 0;
            gas.ElectronDiffusion(0.0, avg_e, 0.0, 0., 0., 0., dl, dt);
            grDiff->SetPoint(graph_idx, i, dt);








            // NEW: Record Drift Velocity at avg_e
            double vx = 0., vy = 0., vz = 0.;
            gas.ElectronVelocity(0.0, avg_e, 0.0, 0., 0., 0., vx, vy, vz);
            double v_mag = std::sqrt(vx*vx + vy*vy + vz*vz);
            grVel->SetPoint(graph_idx, i, v_mag);








            TGraph* grSpatial = new TGraph();
            int pt_idx = 0;
















            
            const int nY = int(std::lround((Y_MAX - Y_MIN) / y_step)) + 1;
            for (int iy = 0; iy < nY; ++iy) {
                const double y_spawn = Y_MIN + y_step * iy;








                const int nElectrons = 100;
                int leaked = 0;








#ifdef _OPENMP
               
                #pragma omp parallel
                {
                    AvalancheMC td;
                    td.SetSensor(&sensor);
                    td.SetDistanceSteps(0.02);
                    td.SetTimeWindow(0.0, 40000.0);




                    #pragma omp for schedule(dynamic, 16) reduction(+:leaked)
                    for (int e = 0; e < nElectrons; ++e) {
                        const double x0 = -1.4 + (2.8 * e) / (nElectrons - 1.0);
                        td.DriftElectron(x0, y_spawn, 0.005, 0.0);




                        const int nPoints = td.GetNumberOfDriftLinePoints();
                        double x1 = 0, y1 = 0, z1 = 0, t1 = 0;
                        if (nPoints > 0) td.GetDriftLinePoint(nPoints - 1, x1, y1, z1, t1);




                        if (nPoints > 0 && (y1 < -15.98) && (std::abs(x1) < 1.5)) {
                            leaked += 1;
                        } else {
                            #pragma omp critical(lostfile)
                            lostFile << y_spawn << "\t" << x0 << "\t" << x1 << "\t"
                                     << y1 << "\t" << t1 << "\t" << nPoints << "\n";
                        }
                    }
                }




#else
                for (int e = 0; e < nElectrons; ++e) {
                    const double x0 = -1.4 + (2.8 * e) / (nElectrons - 1.0);
                    local_drift.DriftElectron(x0, y_spawn, 0.005, 0.0);




                    const int nPoints = local_drift.GetNumberOfDriftLinePoints();
                    double x1 = 0, y1 = 0, z1 = 0, t1 = 0;
                    if (nPoints > 0) local_drift.GetDriftLinePoint(nPoints - 1, x1, y1, z1, t1);




                    if (nPoints > 0 && (y1 < -15.9) && (std::abs(x1) < 1.5)) {
                        leaked++;
                    } else {
                        lostFile << y_spawn << "\t" << x0 << "\t" << x1 << "\t"
                                 << y1 << "\t" << t1 << "\t" << nPoints << "\n";
                    }
                }
#endif












                double pass_pct = 100.0 * leaked / nElectrons;
                outFile << y_spawn << "\t" << pass_pct << "\n";
                grSpatial->SetPoint(pt_idx++, y_spawn, pass_pct);
                std::cout << "  Y = " << y_spawn << " -> " << pass_pct
                          << "% passed\n" << std::flush;
            }








            outFile.close();








            grSpatial->SetLineColor(colors[graph_idx]);
            grSpatial->SetLineWidth(2);
            grSpatial->SetMarkerStyle(20);
            grSpatial->SetMarkerSize(0.25);      // quarter of the ROOT default
            grSpatial->SetMarkerColor(colors[graph_idx]);
            mg->Add(grSpatial, "LP");
            legend->AddEntry(grSpatial, Form("%d%% CF4", i), "lp");








            // INDIVIDUAL FITS & PLOTS
            TGraph* grForFit = (TGraph*)grSpatial->Clone(Form("grFit_%d", i));
            TF1* fitFunc = new TF1(Form("fit_%d", i), "[0]/(x - [1]) + [2]", -15.0, 15.0);
            fitFunc->SetParameters(10.0, -20.0, 50.0);
            fitFunc->SetLineColor(colors[graph_idx]);
            grForFit->Fit(fitFunc, "Q");








            TCanvas* cSingle = new TCanvas(Form("cSingle_%d", i), Form("CF4 %d%% Efficiency", i), 800, 600);
            cSingle->SetGrid();
            grForFit->SetTitle(Form("Spatial Efficiency: %d%% CF4;Spawn Y Position [cm];Electrons Passed [%%]", i));
            grForFit->Draw("AP");
            cSingle->SaveAs((out_dir + "/" + std::to_string(i) + "_EPass_v_Y.pdf").c_str());
            delete cSingle;








            TGraph* grCurveOnly = new TGraph(fitFunc);
            grCurveOnly->SetLineColor(colors[graph_idx]);
            grCurveOnly->SetLineWidth(2);
            mgCurves->Add(grCurveOnly, "L");
           
            graph_idx++;
        }








        // ==========================================================
        // 4. FINAL OUTPUT PLOTS
        // ==========================================================








        // Plot Spatial Efficiency TOTAL
        TCanvas* cEff = new TCanvas("cEff", "Total Spatial Efficiency", 1000, 600);
        cEff->SetRightMargin(0.22);
        cEff->SetGrid();
        mg->SetMinimum(0.0);
        mg->SetMaximum(100.0);
        mg->Draw("A");
        
        mg->GetXaxis()->SetLimits(Y_MIN, Y_MAX);
        mg->GetXaxis()->SetNdivisions(510);
        gPad->Modified();
        gPad->Update();
        legend->Draw();
        cEff->SaveAs((out_dir + "/Total_EPass_v_Y.pdf").c_str());








        // Plot Spatial Efficiency TOTAL CURVES
        TCanvas* cCurves = new TCanvas("cCurves", "Total Fit Curves", 1000, 600);
        cCurves->SetRightMargin(0.22);
        cCurves->SetGrid();
        mgCurves->SetMinimum(0.0);
        mgCurves->SetMaximum(100.0);
        mgCurves->Draw("A");
        mgCurves->GetXaxis()->SetLimits(Y_MIN, Y_MAX);
        mgCurves->GetXaxis()->SetNdivisions(510);
        gPad->Modified();
        gPad->Update();
        legend->Draw();
        cCurves->SaveAs((out_dir + "/TotalCurves_EPass_v_Y.pdf").c_str());








        // Plot Transverse Diffusion
        TCanvas* cDiff = new TCanvas("cDiff", "Transverse Diffusion", 800, 600);
        cDiff->SetGrid();
        grDiff->Draw("APL");
        cDiff->SaveAs((out_dir + "/strip_transverse_diffusion_vs_cf4.pdf").c_str());








        // NEW: Plot Drift Velocity
        TCanvas* cVelSweep = new TCanvas("cVelSweep", "Drift Velocity vs CF4", 800, 600);
        cVelSweep->SetGrid();
        grVel->Draw("APL");
        cVelSweep->SaveAs((out_dir + "/strip_drift_velocity_vs_cf4.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/strip_drift_velocity_vs_cf4.pdf\n";
    } else if (switch_mode == 2) {
        // ==================================================================
        //  MODE 2 - TRANSVERSE / LONGITUDINAL DIFFUSION MEASUREMENT
        //
        //  Launch N electrons from a single point on the axis and follow the
        //  cloud down the cage, sampling the beam width at fixed heights.
        //
        //  The key idea: do NOT just look at where the electrons land. Each
        //  drift line already contains the full trajectory, so one set of N
        //  electrons gives sigma(L) at EVERY depth for free. Endpoint-only
        //  would need a separate run per drift length.
        //
        //  Diffusion theory says   sigma_x = D_T * sqrt(L)
        //  so sigma_x plotted against sqrt(L) is a straight line through the
        //  origin whose SLOPE is D_T in cm^(1/2). That slope, compared with
        //  what MediumGas tabulates, validates the whole transport chain.
        // ==================================================================
        // The gas table is loaded inside mode 0 and mode 1, so mode 2 must
        // load its own. Without this the MediumGas is empty and every
        // transport number comes back zero.
        const std::string gas_file_m2 = "gas_tables/cf4_80_iso_20.gas";
        if (!gas.LoadGasFile(gas_file_m2)) {
            std::cerr << "FATAL: could not load " << gas_file_m2 << "\n";
            gSystem->Exit(1);
        }
        std::cout << "-> Gas: " << gas_file_m2 << "\n";


        const int    nElectrons = 100000;
        const double x_start    =  0.0;    // on axis
        const double y_start    =  0.0;    // exact middle of the cage
        const double z_start    =  0.005;
        const double y_stop     = -15.0;   // deepest sampling plane
        const double dy_check   =  0.5;    // sampling plane spacing [cm]


        // Average |E| over the drift volume, computed here because mode 1's
        // copy lives inside its own block and is not visible from here.
        // Same scan window as mode 1 so the two quote the same field.
        double sum_e2 = 0.0;
        int n_pts2 = 0;
        for (double x = -1.490; x <= 1.491; x += 0.1) {
            for (double y = -15.0; y <= 15.0; y += 0.2) {
                double ex = 0., ey = 0., ez = 0., v = 0.;
                Medium* m = nullptr;
                int status = 0;
                elm.ElectricField(x, y, 0.005, ex, ey, ez, v, m, status);
                if (status == 0) {
                    sum_e2 += std::sqrt(ex * ex + ey * ey + ez * ez);
                    ++n_pts2;
                }
            }
        }
        const double avg_e = (n_pts2 > 0) ? (sum_e2 / n_pts2) : 33.0;
        std::cout << "-> Average field for mode 2: " << avg_e << " V/cm\n";


        const int nCheck = int(std::lround((y_start - y_stop) / dy_check));
        std::vector<double> yCheck(nCheck);
        for (int c = 0; c < nCheck; ++c) yCheck[c] = y_start - dy_check * (c + 1);


        std::vector<long long> nAt(nCheck, 0);
        std::vector<double> sX(nCheck, 0.0), sX2(nCheck, 0.0);
        std::vector<double> sT(nCheck, 0.0), sT2(nCheck, 0.0);
        std::vector<double> sX4(nCheck, 0.0);   // for the tail/kurtosis check


        // Draw the first nPlot trajectories only. ViewDrift stores every
        // step of every line it is given: 1000 lines x ~3000 points would
        // make a PDF nothing can open. Plotting is switched off after
        // nPlot so the remaining electrons still count toward the
        // statistics without being drawn.
        const int nPlot = 100;


        AvalancheMC dif;
        dif.SetSensor(&sensor);
        dif.SetDistanceSteps(0.01);
        dif.SetTimeWindow(0.0, 200000.0);   // deliberately generous: a
                                            // truncated line would silently
                                            // bias sigma at large L
        dif.EnablePlotting(&driftView);


        TH1D* hEnd = new TH1D("hEnd", "Final x after full drift;x [cm];Electrons",
                              150, -1.5, 1.5);


        std::cout << "\nMODE 2: " << nElectrons << " electrons from ("
                  << x_start << ", " << y_start << ")\n"
                  << "Sampling " << nCheck << " planes every " << dy_check << " cm\n";


        int nComplete = 0, nAbsorbed = 0;
        for (int e = 0; e < nElectrons; ++e) {
            if (e == nPlot) dif.DisablePlotting();
            dif.DriftElectron(x_start, y_start, z_start, 0.0);
            const int nPts = dif.GetNumberOfDriftLinePoints();
            if (nPts < 2) continue;


            double xE, yE, zE, tE;
            dif.GetDriftLinePoint(nPts - 1, xE, yE, zE, tE);
            if (yE < y_stop) { ++nComplete; hEnd->Fill(xE); }
            else             { ++nAbsorbed; }


            // Walk the trajectory once, crossing sampling planes in order.
            int c = 0;
            double xP, yP, zP, tP;
            dif.GetDriftLinePoint(0, xP, yP, zP, tP);
            for (int p = 1; p < nPts && c < nCheck; ++p) {
                double x, y, z, t;
                dif.GetDriftLinePoint(p, x, y, z, t);
                while (c < nCheck && yP > yCheck[c] && y <= yCheck[c]) {
                    const double dyseg = yP - y;
                    const double f = (dyseg > 1e-12) ? (yP - yCheck[c]) / dyseg : 0.0;
                    const double xi = xP + f * (x - xP);
                    const double ti = tP + f * (t - tP);
                    nAt[c] += 1;
                    sX[c]  += xi;   sX2[c] += xi * xi;   sX4[c] += xi * xi * xi * xi;
                    sT[c]  += ti;   sT2[c] += ti * ti;
                    ++c;
                }
                xP = x; yP = y; zP = z; tP = t;
            }
            if ((e + 1) % 100 == 0)
                std::cout << "  " << (e + 1) << " / " << nElectrons << "\r" << std::flush;
        }
        std::cout << "\nReached the bottom: " << nComplete
                  << "   absorbed early: " << nAbsorbed << "\n";


        // ---- gas prediction for comparison ----
        double dl_tab = 0.0, dt_tab = 0.0;
        gas.ElectronDiffusion(0.0, avg_e, 0.0, 0., 0., 0., dl_tab, dt_tab);
        double vx_t = 0., vy_t = 0., vz_t = 0.;
        gas.ElectronVelocity(0.0, avg_e, 0.0, 0., 0., 0., vx_t, vy_t, vz_t);
        const double v_tab = std::sqrt(vx_t * vx_t + vy_t * vy_t + vz_t * vz_t);


        std::ofstream df(out_dir + "/diffusion_profile.txt");
        df << "L_cm\tn\tmean_x\tsigma_x\tDT_local\tmean_t_ns\tsigma_t_ns\tDL_local\texcess_kurt\n";


        TGraph* grSig  = new TGraph();   // sigma_x vs sqrt(L)  -> slope is D_T
        TGraph* grDTL  = new TGraph();   // D_T vs L            -> should be flat
        TGraph* grSigT = new TGraph();   // sigma_t vs sqrt(L)
        int k = 0;


        for (int c = 0; c < nCheck; ++c) {
            if (nAt[c] < 50) continue;
            const double N  = double(nAt[c]);
            const double L  = y_start - yCheck[c];
            const double mx = sX[c] / N;
            double vx = sX2[c] / N - mx * mx;            if (vx < 0) vx = 0;
            const double sx = std::sqrt(vx);
            const double mt = sT[c] / N;
            double vt = sT2[c] / N - mt * mt;            if (vt < 0) vt = 0;
            const double st = std::sqrt(vt);


            // excess kurtosis about 0; 0 = Gaussian, <0 = truncated by the walls
            const double kurt = (vx > 0) ? (sX4[c] / N) / (vx * vx) - 3.0 : 0.0;


            // longitudinal spread converted from time to length
            const double sL   = st * v_tab;
            const double DT_l = sx / std::sqrt(L);
            const double DL_l = sL / std::sqrt(L);


            df << L << "\t" << nAt[c] << "\t" << mx << "\t" << sx << "\t" << DT_l
               << "\t" << mt << "\t" << st << "\t" << DL_l << "\t" << kurt << "\n";


            grSig ->SetPoint(k, std::sqrt(L), sx);
            grSigT->SetPoint(k, std::sqrt(L), sL);
            grDTL ->SetPoint(k, L, DT_l);
            ++k;
        }
        df.close();


        // ---- fit sigma_x = D_T * sqrt(L), forced through the origin ----
        TF1* fT = new TF1("fT", "[0]*x", 0.0, std::sqrt(y_start - y_stop) * 1.05);
        fT->SetParameter(0, dt_tab);
        grSig->Fit(fT, "Q");
        const double DT_fit = fT->GetParameter(0);
        const double DT_err = fT->GetParError(0);


        TF1* fL = new TF1("fL", "[0]*x", 0.0, std::sqrt(y_start - y_stop) * 1.05);
        fL->SetParameter(0, dl_tab);
        grSigT->Fit(fL, "Q");
        const double DL_fit = fL->GetParameter(0);


        std::cout << "\n===== DIFFUSION SUMMARY (E = " << avg_e << " V/cm) =====\n";
        std::cout << "  measured D_T = " << DT_fit << " +/- " << DT_err << " cm^1/2\n";
        std::cout << "  tabulated D_T = " << dt_tab << " cm^1/2   ratio "
                  << (dt_tab > 0 ? DT_fit / dt_tab : 0) << "\n";
        std::cout << "  measured D_L = " << DL_fit << " cm^1/2\n";
        std::cout << "  tabulated D_L = " << dl_tab << " cm^1/2   ratio "
                  << (dl_tab > 0 ? DL_fit / dl_tab : 0) << "\n";
        std::cout << "  sigma_x over the measured " << (y_start - y_stop) << " cm = "
                  << DT_fit * std::sqrt(y_start - y_stop) << " cm\n";
        std::cout << "  extrapolated to a full 30.9 cm drift = "
                  << DT_fit * std::sqrt(30.9) << " cm\n";
        std::cout << "  drift velocity (table) = " << v_tab << " cm/ns\n";


        // ---- plots ----
        TCanvas* cDr = new TCanvas("cDr", "Diffusion Drift Lines", 600, 600);
        driftView.SetCanvas(cDr);
        driftView.Plot(true);
        TLine* dL = new TLine(-1.5, -16.0, -1.5, 16.0);
        TLine* dR = new TLine( 1.5, -16.0,  1.5, 16.0);
        dL->SetLineColor(kGray + 2); dL->SetLineWidth(2); dL->Draw();
        dR->SetLineColor(kGray + 2); dR->SetLineWidth(2); dR->Draw();
        cDr->SaveAs((out_dir + "/diffusion_drift_lines.pdf").c_str());
        std::cout << "-> Saved " << out_dir << "/diffusion_drift_lines.pdf ("
                  << nPlot << " of " << nElectrons << " lines drawn)\n";


        TCanvas* cH = new TCanvas("cH", "Final x", 800, 600);
        cH->SetGrid();
        hEnd->SetLineColor(kBlue + 1); hEnd->SetLineWidth(2);
        hEnd->Fit("gaus", "Q");
        hEnd->Draw();
        cH->SaveAs((out_dir + "/diffusion_endpoint_x.pdf").c_str());


        TCanvas* cS = new TCanvas("cS", "sigma vs sqrt(L)", 800, 600);
        cS->SetGrid();
        grSig->SetTitle("Transverse spread;#sqrt{L} [cm^{1/2}];#sigma_{x} [cm]");
        grSig->SetMarkerStyle(20); grSig->SetMarkerSize(0.5);
        grSig->SetMarkerColor(kRed + 1); grSig->SetLineColor(kRed + 1);
        grSig->Draw("AP");
        cS->SaveAs((out_dir + "/diffusion_sigma_vs_sqrtL.pdf").c_str());


        TCanvas* cD = new TCanvas("cD", "D_T vs L", 800, 600);
        cD->SetGrid();
        grDTL->SetTitle("Local D_{T} (flat = pure diffusion);Drift length L [cm];#sigma_{x}/#sqrt{L} [cm^{1/2}]");
        grDTL->SetMarkerStyle(20); grDTL->SetMarkerSize(0.5);
        grDTL->SetMarkerColor(kGreen + 2); grDTL->SetLineColor(kGreen + 2);
        grDTL->Draw("APL");
        TLine* lTab = new TLine(0.0, dt_tab, y_start - y_stop, dt_tab);
        lTab->SetLineColor(kBlack); lTab->SetLineStyle(2); lTab->Draw();
        cD->SaveAs((out_dir + "/diffusion_DT_vs_L.pdf").c_str());


        std::cout << "-> Wrote diffusion_profile.txt and three PDFs to " << out_dir << "\n";
    }








    gSystem->Exit(0);
    return 0;
}


