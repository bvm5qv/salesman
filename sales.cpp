#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TObject.h"

typedef struct {
  double lon, lat; //in degrees
} COORD;

constexpr int NMAX = 2500;          //Max amount of cities allowed
constexpr double g_R = 6371;        //Earth's radius in km
constexpr double g_Pi = 3.14159265; //pi

//project degrees latitude to mercator latitude
double MercatorY(double lat_deg)
{
  double lat = lat_deg * g_Pi/180;
  return TMath::Log(TMath::Tan(g_Pi/4 + lat/2));
}

//make Mercator projection plot
TH2F *MercatorMap(double latMin, double latMax)
{
  TH2F *hm = new TH2F("hm", "Mercator", 180, -180, 180,
                      160, MercatorY(latMin), MercatorY(latMax));

  std::ifstream in("earth.dat");
  if(!in.is_open())
  {
    Error("MercatorMap", "Cannot open earth.dat");
    return nullptr;
  }

  double lon, lat;
  while(in >> lon >> lat)
  {
    if(lat < latMin || lat > latMax) continue;
    hm->Fill(lon, MercatorY(lat));
  }

  in.close();
  return hm;
}

//draw info box on final plot
void DrawInfoBox(TCanvas *c, double initPathLen, double finalPathLen, double compTime)
{
  c->cd();
  TLegend *leg = new TLegend(0.65, 0.65, 0.9, 0.85);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetTextSize(0.03);
  leg->SetMargin(0.1);
  leg->SetTextFont(42);

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(2);
  ss << "Initial Path Length: " << initPathLen << "km\n";
  leg->AddEntry((TObject*)0, ss.str().c_str(), "");

  ss.str(""); ss.clear();
  ss << "Final Path Length: " << finalPathLen << "km\n";
  leg->AddEntry((TObject*)0, ss.str().c_str(), "");

  ss.str(""); ss.clear();
  ss << "Runtime: " << compTime << "s\n";
  leg->AddEntry((TObject*)0, ss.str().c_str(), "");

  leg->Draw();
  c->Update();
}

//fill the array of city locations
int GetData(char* fname, COORD *cities){
  FILE* fp=fopen(fname,"r");
  const int bufsiz=1000;
  char line[bufsiz+1];
  int ncity=0;
  while(1){
    fgets(line,bufsiz,fp);
    if (line[0]=='#') continue;  // skip comments
    if (feof(fp)) break;
    // we only scan for two numbers at start of each line
    sscanf(line,"%lf %lf",&cities[ncity].lon,&cities[ncity].lat);    
    ncity++;
  }
  fclose(fp);
  return ncity;
}

//find random double between 0 and 1
double GetRandProb()
{ 
  static std::random_device rd;   //non-deterministic seed
  static std::mt19937 gen(rd());  //Mersenne Twister RNG
  static std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(gen); 
}

//find random int between 0 and n (inclusive)
int GetRandInt(int n)
{
  static std::random_device rd;   //non-deterministic seed
  static std::mt19937 gen(rd());  //Mersenne Twister RNG
  std::uniform_int_distribution<> dist(0,n);
  return dist(gen);
}

//calculate distance between two points on Earth
double GetDist(const COORD pos1, const COORD pos2)
{
  double dlat = pos2.lat - pos1.lat;
  double dlon = pos2.lon - pos1.lon;
  double a = std::sin(dlat/2 * g_Pi/180)*std::sin(dlat/2 * g_Pi/180) + std::cos(pos1.lat * g_Pi/180) 
           * std::cos(pos2.lat * g_Pi/180) * std::sin(dlon/2 * g_Pi/180)*std::sin(dlon/2 * g_Pi/180);
  double c = 2*std::atan2(std::sqrt(a), std::sqrt(1 - a));
  return g_R*c;
}

//calculate path length of full circuit
double PathLength(const COORD *cities, const int nCities)
{
  double length = 0;
  for(int i = 0; i < nCities - 1; i++)
    length += GetDist(cities[i], cities[i+1]);
  length += GetDist(cities[nCities - 1], cities[0]);
  return length;
} 

//returns change in path length from permutation of two cities
double TwoOptDiff(const COORD *cities, const int nCities, const int i, const int j)
{
  int in = (i+1==nCities) ? 0 : i+1;
  int jn = (j+1==nCities) ? 0 : j+1;
  if(in == j) return 0.0;
  if(i == 0 && j == nCities - 1) return 0.0;

  double len0 = GetDist(cities[i], cities[in]) + GetDist(cities[j], cities[jn]);
  double len1 = GetDist(cities[i], cities[j]) + GetDist(cities[in], cities[jn]);
  return len1 - len0;
}

void ApplyTwoOpt(COORD *cities, const int i, const int j)
{
  int a = i + 1;
  int b = j;
  while(a < b)
    std::swap(cities[a++], cities[b--]);
}

//Perform simulated annealing
void Anneal(COORD *cities, const int nCities, const double T0, const double Tmin, const double alpha, std::vector<double> &vTemp, std::vector<double> &vLen) 
{
  double currentLength = PathLength(cities, nCities);
  double bestLength = currentLength;
  COORD newPath[NMAX];
  for(int i = 0; i < nCities; i++) newPath[i] = cities[i];

  double T = T0;
  vTemp.push_back(T);
  vLen.push_back(currentLength);
  while(T > Tmin)
  {
    for(int i = 0; i < 10*nCities; i++)
    {
      int j = GetRandInt(nCities-3);
      int k = GetRandInt(nCities-1);
      if(j > k) std::swap(j, k);
      if(k <= j + 1) continue;
      if(j == 0 && k == nCities - 1) continue;

      double dLen = TwoOptDiff(newPath, nCities, j, k);
      if(dLen <= 0.0 || GetRandProb() < std::exp(-dLen/T))
      {
        ApplyTwoOpt(newPath, j, k);
        currentLength += dLen;
        if(currentLength < bestLength)
        {
          bestLength = currentLength;
          for(int l = 0; l < nCities; l++) cities[l] = newPath[l];
        }
      }
    }
    T *= alpha; //update temperature
    vTemp.push_back(T);
    vLen.push_back(currentLength);
  }
}

//========================================================================================================================
//========================================================================================================================
//========================================================================================================================
int main(int argc, char *argv[])
{
  COORD cities[NMAX];
  COORD earth[21450];
  if (argc<5){
    printf("Please provide a data file path, initial temperature, minimum temperature, and cooling constant (between 0 and 1) as arguments\n");
    return 1;
  }

  const auto t0 = std::chrono::system_clock::now();
  int nCities=GetData(argv[1],cities);
  printf("Read %d cities from data file\n",nCities);

  double T0 = std::stod(argv[2]);
  double Tmin = std::stod(argv[3]);
  double alpha = std::stod(argv[4]);
  double initLen = PathLength(cities, nCities);

  std::vector<double> vTemp{};
  std::vector<double> vLen{};
  Anneal(cities, nCities, T0, Tmin, alpha, vTemp, vLen);
  
  double finalLen = PathLength(cities, nCities);
  const auto t1 = std::chrono::system_clock::now();
  const auto dt = std::chrono::duration<double>(t1 - t0);
  std::cout << "Finished computation in " << dt.count() << " seconds.\n";
  std::cout << "Initial Path Length: " << initLen << " km\nFinal Path Length: " << finalLen << " km\n";

//========================================================================================================================
//=============================================== Graphing ===============================================================
//========================================================================================================================

  TApplication App("App", &argc, argv);

  //plot annealing schedule
  TGraph *tgA = new TGraph(vTemp.size(), vTemp.data(), vLen.data());
  tgA->SetTitle("Annealing Schedule;Temperature(km);Path Length(km)");
  tgA->SetMarkerStyle(20);
  tgA->SetMarkerColor(kBlue);
  tgA->SetLineColor(kBlue);
  TCanvas *c1 = new TCanvas("c1", "canvas", 800, 600);
  tgA->Draw("ALP");

  //plot path
  double lon[nCities+1];
  double lat[nCities+1];
  for(int i = 0; i < nCities; i++)
  {
    lon[i] = cities[i].lon;
    lat[i] = MercatorY(cities[i].lat);
  }
  lon[nCities] = lon[0];
  lat[nCities] = lat[0];
  TGraph *tgP = new TGraph(nCities+1, lon, lat);
  tgP->SetTitle("Final Path;longitude;latitude");
  tgP->SetMarkerStyle(2);
  tgP->SetMarkerColor(kRed);
  tgP->SetLineColor(kRed);
  tgP->SetLineWidth(2);
  TH2F *hm = MercatorMap(-80.5, 80.5);
  TCanvas *c2 = new TCanvas("c2", "canvas", 800, 600);
  hm->Draw("COLZ");
  tgP->Draw("L SAME");
  DrawInfoBox(c2, initLen, finalLen, dt.count());

//========================================================================================================================
//=============================================== Output File ============================================================
//========================================================================================================================

  if(argc >= 5)
  {
    std::ofstream fout;
    fout.open(argv[5], std::ios::out | std::ios::trunc);
    assert(fout.is_open() && "error opening output file");
    for(int i = 0; i < nCities; i++)
      fout << cities[i].lon << ' ' << cities[i].lat << '\n'; 
    fout.close();
  }

  std::cout << "Press ^c to exit\n";
  App.SetIdleTimer(1800,".q");  // set up a failsafe timer to end the program  
  App.Run();

  return 0;
}
