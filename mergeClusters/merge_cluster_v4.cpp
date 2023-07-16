// mergeClusters.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#define MAX_TOURS 100
#define SENSOR 300

struct point
{
public:
	int x, y;
	point(int a = 0, int b = 0)
	{
		x = a;
		y = b;
	}
	point operator+(point second)
	{
		point temp;
		temp.x = x + second.x;
		temp.y = y + second.y;
		return temp;
	}
	int operator==(point second)
	{
		return (second.x == x && second.y == y);
	}
	int distance(point dst)
	{
		return (abs(x - dst.x) + abs(y - dst.y));
	}
	bool IsOnPath(point dst, point vec)
	{
		point next = *this + vec;
		if (distance(dst) > next.distance(dst))
			return true;
		return false;
	}
};

int dis(point a, point b)
{
	int dis = fabs(a.x - b.x) + fabs(a.y - b.y);
	return dis;
}

float cost(vector<point> tour1, vector<point> tour2, int it1, int it2, bool bCross)
{
	float res;
	//    bCross = false                      bCross = true
	//    |-----------|              	//    |-------|
	//    |           |				 	//    |       |
	//   pt1B---------pt1E			 	//   pt1B----pt1E
	//    |           |				 	//       \  /       
	//    |           |				 	//       /  \        
	//   pt2B---------pt2E			 	//   pt2B----pt2E
	//    |           |				 	//    |       |
	//    |-----------|				 	//    |-------|
	point pt1B = tour1[it1];
	point pt1E = tour1[(it1 + 1) % tour1.size()];
	point pt2B = tour2[it2];
	point pt2E = tour2[(it2 + 1) % tour2.size()];

	if (!bCross)
		res = dis(pt1B, pt2B) + dis(pt1E, pt2E) - dis(pt1B, pt1E) - dis(pt2B, pt2E);
	else
		res = dis(pt1B, pt2E) + dis(pt1E, pt2B) - dis(pt1B, pt1E) - dis(pt2B, pt2E);
	return res;
}

struct pivot
{
	int n1, n2, n3;
	bool bCross;
	float cost;
public:
	pivot() {};
	pivot(int i, int j, bool b, float c)
	{
		n1 = i;
		n2 = j;
		bCross = b;
		cost = c;
	}
	pivot(int i, int j, int k, bool b, float c)
	{
		n1 = i;
		n2 = j;
		n3 = k;
		bCross = b;
		cost = c;
	}
};

pivot FindMinCycle(vector<point> tour1, vector<point> tour2)
{
	pivot minP;
	minP.cost = 1000000;
	for (int i = 0; i < tour1.size(); i++)
		for (int j = 0; j < tour2.size(); j++)
		{
			float f = cost(tour1, tour2, i, j, false);
			if (f < minP.cost)
				minP = pivot(i, j, false, f);
			f = cost(tour1, tour2, i, j, true);
			if (f < minP.cost)
				minP = pivot(i, j, true, f);
		}
	return minP;
}

vector<point> merge(vector<point> tour1, vector<point> tour2, pivot p)
{
	//    bCross = false                      bCross = true
	//    |-----------|              	//    |-------|
	//    |           |				 	//    |       |
	//   pt1B---------pt1E			 	//   pt1B----pt1E
	//    |           |				 	//       \  /       
	//    |           |				 	//       /  \        
	//   pt2B---------pt2E			 	//   pt2B-----pt2E
	//    |           |				 	//    |        |
	//    |-----------|				 	//    |--------|
	vector<point> res;
	int it1 = p.n1;
	int it2 = p.n2;
	point pt1B = tour1[it1];
	point pt1E = tour1[(it1 + 1) % tour1.size()];
	point pt2B = tour2[it2];
	point pt2E = tour2[(it2 + 1) % tour2.size()];
	if (p.bCross == false)
	{
		res.push_back(pt1B);
		res.push_back(pt2B);
		for (int i = 1; i <= tour2.size() - 2; i++)
		{
			int j = (it2 - i + tour2.size()) % tour2.size();
			res.push_back(tour2[j]);
		}
		res.push_back(pt2E);
		res.push_back(pt1E);
		for (int i = tour1.size() - 2; i >= 1; i--)
		{
			int j = (it1 - i + tour1.size()) % tour1.size();
			res.push_back(tour1[j]);
		}
	}
	else
	{
		res.push_back(pt1B);
		res.push_back(pt2E);
		for (int i = 1; i <= tour2.size() - 2; i++)
		{
			int j = (it2 + 1 + i) % tour2.size();
			res.push_back(tour2[j]);
		}
		res.push_back(pt2B);
		res.push_back(pt1E);
		for (int i = 1; i <= tour1.size() - 2; i++)
		{
			int j = (it1 - i + tour1.size()) % tour1.size();
			res.push_back(tour1[j]);
		}
	}
	return res;
}

int findSensor(point a, point b[SENSOR])
{
	int index = 9999;
	for (int i = 0; i < SENSOR; i++)
	{
		if (a == b[i]) index = i;
	}
	return index;
}
void  RemoveCross(vector<point>* tour)
{
	if (tour->size() < 4)
		return;
	for (int i = 0; i < tour->size() - 1; i++)
	{
		point p[4];
		for (int j = 0; j < 4; j++)
			p[j] = tour->at((i + j) % tour->size());
		int diff = dis(p[0], p[1]) + dis(p[2], p[3]) - dis(p[0], p[2]) - dis(p[1], p[3]);
		if (diff > 0)
		{
			// change node p[i+1] and p[i+2]
		std:iter_swap(tour->begin() + (i + 1) % tour->size(), tour->begin() + (i + 2) % tour->size());
			printf("\nChanging cross edges for tour at pos %d profit = %d\n", i, diff);
		}
	}
}
void  RemoveLoop(vector<point>* tour)
{
	for (int i = 0; i < tour->size() - 1; i++)
	{
		if (tour->at(i) == tour->at(i + 1))
			tour->erase(tour->begin() + i);
	}
	// for begin and end
	if (tour->at(0) == tour->at(tour->size() - 1))
		tour->pop_back();
}
bool  RemoveCross2(vector<point>* tour)
{
	if (tour->size() < 4)
		return false;
	pivot piv;
	piv.cost = 0;
	for (int i = 0; i < tour->size(); i++)
	{
		point p[2] = { tour->at(i % tour->size()) ,  tour->at((i + 1) % tour->size()) };
		for (int j = i + 1; j < tour->size(); j++)
		{
			point q[2] = { tour->at(j % tour->size()) ,  tour->at((j + 1) % tour->size()) };

			int diff = dis(p[0], p[1]) + dis(q[0], q[1]) - dis(p[0], q[0]) - dis(p[1], q[1]);
			if (diff > piv.cost)
			{
				//std:iter_swap(tour->begin() + (i + 1) % tour->size(), tour->begin() + (i + 2) % tour->size());
				//printf("\n Find cross edges for tour at e(%d,%d)-(%d,%d) with e(%d,%d)-(%d,%d) profit = %d\n",
				//p[0].x, p[0].y, p[1].x, p[1].y, q[0].x, q[0].y, q[1].x, q[1].y, diff);
				piv.n1 = i;
				piv.n2 = j;
				piv.cost = diff;
			}
		}
	}

	if (piv.cost > 0)
	{
		point p[2] = { tour->at(piv.n1 % tour->size()) ,  tour->at((piv.n1 + 1) % tour->size()) };
		point q[2] = { tour->at(piv.n2 % tour->size()) ,  tour->at((piv.n2 + 1) % tour->size()) };
		//printf("\n Apply cross edges(%d, %d) for tour at e(%d,%d)-(%d,%d) with e(%d,%d)-(%d,%d) profit = %2.2f\n",
		//piv.n1, piv.n2, p[0].x, p[0].y, p[1].x, p[1].y, q[0].x, q[0].y, q[1].x, q[1].y, piv.cost);
		vector<point> t1, t2;
		for (int k = 0; k <= piv.n1; k++)
			t1.push_back(tour->at(k));
		for (int k = piv.n1 + 1; k <= piv.n2; k++)
			t2.push_back(tour->at(k));
		for (int k = piv.n2 + 1; k < tour->size(); k++)
			t1.push_back(tour->at(k));
		piv.bCross = false;
		piv.n2 = t2.size() - 1;
		vector<point> res = merge(t1, t2, piv);
		tour->assign(res.begin(), res.end());
		RemoveLoop(tour);
		return true;
	}
	return false;
}



bool  RemoveCross3(vector<point>* tour)
{
	if (tour->size() < 4)
		return false;
	pivot piv;
	piv.cost = 0;

	for (int i = 0; i < tour->size(); i++)
	{
		point p[2] = { tour->at(i % tour->size()) ,  tour->at((i + 1) % tour->size()) };
		for (int j = i + 1; j < tour->size(); j++)
		{
			point q[2] = { tour->at(j % tour->size()) ,  tour->at((j + 1) % tour->size()) };

			for (int k = j + 1; k < tour->size(); k++)
			{
				point r[2] = { tour->at(k % tour->size()) ,  tour->at((k + 1) % tour->size()) };

				int diff0 = dis(p[0], p[1]) + dis(q[0], q[1]) + dis(r[0], r[1]);
				int diff3 = dis(p[0], q[1]) + dis(r[0], p[1]) + dis(q[0], r[1]);

				if (diff0 - diff3 > piv.cost)
				{
					//std:iter_swap(tour->begin() + (i + 1) % tour->size(), tour->begin() + (i + 2) % tour->size());
					//printf("\n Find cross edges for tour at %d, %d, %d with  profit = %d\n",
					//i, j, k, diff0 - diff3);
					piv.n1 = i;
					piv.n2 = j;
					piv.n3 = k;
					piv.cost = diff0 - diff3;
				}
			}
		}
	}

	if (piv.cost > 0)
	{
		point p[2] = { tour->at(piv.n1 % tour->size()) ,  tour->at((piv.n1 + 1) % tour->size()) };
		point q[2] = { tour->at(piv.n2 % tour->size()) ,  tour->at((piv.n2 + 1) % tour->size()) };
		point r[2] = { tour->at(piv.n3 % tour->size()) ,  tour->at((piv.n3 + 1) % tour->size()) };
		//printf("\n Apply cross 3 edges(%d, %d) for tour at e(%d,%d)-(%d,%d) with e(%d,%d)-(%d,%d) and e(%d,%d)-(%d,%d) profit = %2.2f\n",
		//piv.n1, piv.n2, p[0].x, p[0].y, p[1].x, p[1].y,
		//q[0].x, q[0].y, q[1].x, q[1].y,
		//r[0].x, r[0].y, r[1].x, r[1].y,
		//piv.cost);

		vector<point>::const_iterator f0 = tour->begin();
		vector<point>::const_iterator f1 = tour->begin() + piv.n1 + 1;
		vector<point>::const_iterator f2 = tour->begin() + piv.n2 + 1;
		vector<point>::const_iterator f3 = tour->begin() + piv.n3 + 1;
		vector<point>::const_iterator f4 = tour->end();
		vector<point> newVec1(f0, f1);
		vector<point> newVec2(f1, f2);
		vector<point> newVec3(f2, f3);
		vector<point> newVec4(f3, f4);
		newVec1.insert(newVec1.end(), newVec3.begin(), newVec3.end());
		newVec1.insert(newVec1.end(), newVec2.begin(), newVec2.end());
		newVec1.insert(newVec1.end(), newVec4.begin(), newVec4.end());
		tour->clear();
		tour->insert(tour->end(), newVec1.begin(), newVec1.end());
		return true;
	}
	return false;
}




int  Cost(vector<point>* tour)
{
	int res = 0;
	for (int i = 0; i < tour->size(); i++)
	{
		res += dis(tour->at(i), tour->at((i + 1) % tour->size()));
	}
	return res;
}

struct Rec
{
	point p;
	double sumAOI;
	int cnt;
};
float avgAOI(vector<point>* tour)
{
	vector<Rec> res;
	for (int i = 0; i < tour->size(); i++)
	{
		Rec r;
		r.p = tour->at(i);
		r.sumAOI = 0;
		r.cnt = 1;
		bool bfound = false;
		for (int j = 0; j < res.size(); j++)
		{
			if (res.at(j).p == r.p)
			{
				bfound = true;
				res.at(j).cnt++;
				break;
			}
		}
		if (bfound)
			continue;
		int d = 0;
		for (int j = 0; j < tour->size(); j++)
		{
			d += dis(tour->at((i + j) % tour->size()), tour->at((i + j + 1) % tour->size()));
			if (tour->at((i + j + 1) % tour->size()) == r.p)
			{
				r.sumAOI += ((double)d*d - d) / (double)2;
				d = 0;
			}
		}
		res.push_back(r);
	}
	int dTour = Cost(tour);
	float avgAoI = 0;
	for (int i = 0; i < res.size(); i++)
	{
		avgAoI += res.at(i).sumAOI;
	}
	avgAoI /= (dTour * res.size());
	return avgAoI;
}

bool  InsertLoop(vector<point>* tour)
{
	float currentAoI = avgAOI(tour);
	float maxdiff = 0;
	int n1, n2;
	for (int i = 0; i < tour->size(); i++)
	{
		for (int j = 0; j < tour->size(); j++)
		{
			vector<point> t2;
			if (tour->at(j) == tour->at(i) || tour->at(j) == tour->at((i + 1) % tour->size()))
				continue;
			int d1 = tour->at(j).distance(tour->at(i));
			int d2 = tour->at(j).distance(tour->at((i + 1) % tour->size()));
			int d3 = tour->at(i).distance(tour->at((i + 1) % tour->size()));
			if (d1 + d2 > d3 + 20)
				continue;
			t2.insert(t2.end(), tour->begin(), tour->begin() + i + 1);
			t2.insert(t2.end(), tour->at(j));
			t2.insert(t2.end(), tour->begin() + i + 1, tour->end());
			float newAoI = avgAOI(&t2);
			if (currentAoI - newAoI > maxdiff)
			{
				maxdiff = currentAoI - newAoI;
				n1 = i;
				n2 = j;
				//printf("Check : Insertiion of (%d,%d) between %d and %d gain=%2.2f\n",
				//	tour->at(j).x, tour->at(j).y, i, i + 1, maxdiff);
				//printf("d1 = %d d2= %d d3=%d\n", d1, d2, d3);
			}
		}
	}
	if (maxdiff > 0)
	{
		vector<point> t2;
		t2.insert(t2.end(), tour->begin(), tour->begin() + n1 + 1);
		t2.insert(t2.end(), tour->at(n2));
		t2.insert(t2.end(), tour->begin() + n1 + 1, tour->end());
		printf("Inserting (%d,%d) between %d and %d \n",
			tour->at(n2).x, tour->at(n2).y, n1, n1 + 1);
		tour->clear();
		tour->insert(tour->end(), t2.begin(), t2.end());
		return true;
	}
	return false;
}

int main()
{

	ofstream costs;

	int sizeOfClusters[MAX_TOURS];
	int nSensor, nCluster, m, nUAV;
	cout << "Number of Sensors : \n";
	cin >> nSensor;
	cout << "Number of clusters : \n";
	cin >> nCluster;
	cout << "Method ? averageAOI = 1 , hopCount = 2 : \n";
	cin >> m;
	cout << "Number of UAV : \n";
	cin >> nUAV;

	int m1 = -1;
	if (nCluster == 1)
	{
		cout << "proposed or opt ? proposed = 1 , opt = 2 : \n";
		cin >> m1;
	}

	string method;
	if (m == 1) method = "AOISum"; else method = "HopCount";

	string dists = "costs - " + to_string(nSensor) + " - " + to_string(nUAV) + "UAV.csv";
	costs.open(dists);
	for (int i = 0; i < nUAV; i++) costs << "Dist UAV" << i + 1 << ", " << "AvgAOI UAV" << i + 1 << ", " << "nCluster UAV" << i + 1 << ", " << "nSensor UAV" << i + 1 << ", ";
	costs << endl;

	int nCluster3 = nCluster - 1;
	while (nCluster3 < 20)
	{

		int sensorIndexs[500];
		for (int i = 0; i < 500; i++) sensorIndexs[i] = -1;
		nCluster = ++nCluster3;
		cout << " ====> nCluster = " << nCluster << endl << endl;

		string fileName, sensorFileName, sensorLocation;
		switch (nSensor)
		{
		case 10:
			fileName = "test10";
			break;
		case 51:
			fileName = "eil51";
			break;
		case 76:
			fileName = "eil76";
			break;
		case 101:
			fileName = "eil101";
			break;
		case 52:
			fileName = "berlin52";
			break;
		case 48:
			fileName = "att48";
			break;
		case 280:
			fileName = "a280";
			break;
		default:
			break;
		}
		sensorFileName = fileName + ".tsp";
		sensorLocation = fileName + ".csv";



		vector<point> tours[MAX_TOURS];

		string s[MAX_TOURS];

		if (nCluster == 1)
		{
			if (m1 == 1) s[0] = sensorFileName + " - " + method + ".opt.tour";
			else s[0] = sensorFileName + ".opt.tour";
		}
		else {

			for (int i = 0; i < nCluster; i++)
			{
				s[i] = to_string(nSensor) + "Sensors - " + to_string(nCluster) + "cluster - " + to_string(i + 1) + ".tsp - " + method + ".opt.tour";
			}
		}
		int index = 0, index1;
		ifstream tour;
		string temp;
		for (int i = 0; i < nCluster; i++)
		{
			tour.open(s[i]);
			for (int i = 0; i < 5; i++) getline(tour, temp);
			int j = 0;
			tour >> index; tour >> index1;
			while (index != -1)
			{
				tours[i].push_back(point(index, index1));
				tour >> index; tour >> index1;
				j++;
			}
			tour.close();
		}

		/**
		tours[0] = { point(415, 635), point(345, 750), point(145, 665), point(300, 465), point(420, 555), point(565, 575), point(575, 665), point(595, 360), point(660, 180), point(410, 250), point(480, 415), point(560, 365), point(520, 585) };
		tours[1] = { point(1170, 65), point(1215, 245), point(1320, 315), point(1465, 200), point(1530, 5), point(1740, 245) };
		tours[2] = { point(1250, 400), point(1220, 580), point(1340, 725), point(1605, 620) };
		tours[3] = { point(555, 815), point(510, 875), point(475, 960), point(525, 1000), point(580, 1175), point(650, 1130), point(875, 920) };
		tours[4] = { point(725, 370), point(700, 500), point(700, 580), point(685, 595), point(685, 610), point(720, 635), point(760, 650), point(795, 645), point(845, 655), point(845, 680), point(880, 660), point(845, 655), point(835, 625), point(830, 610), point(770, 610), point(685, 610), point(605, 625), point(685, 610), point(685, 595), point(700, 580), point(700, 500), point(830, 485), point(975, 580), point(945, 685) };
		tours[5] = { point(25, 185), point(25, 230), point(95, 260) };
		tours[6] = { point(1150, 1160) };
		*/

		// ToDo : copy tour1,2,3,4 => tours
		int nTour1 = nCluster;
		int nTour = nCluster;
		//merging the tours until becomes one tour.

		ofstream plotSensor;
		string sensorsPlot = to_string(nSensor);
		sensorsPlot += "sensors - " + to_string(nCluster) + "cluster.csv";
		plotSensor.open(sensorsPlot);
		for (int i = 0; i < nTour; i++)
		{
			if (tours[i].size() == 1)
			{
				tours[i].push_back(tours[i].front());
			}
		}
		int level = 0;
		int ns = 0;
		int final_minP = 0;
		int repeat = nUAV;
		int tourNumber = 1;
		int jTour = 0;
		int boundaryTour = 0;
		int visited[100];
		for (int i = 0; i < nCluster; i++) visited[i] = -1;

		while (repeat > 0)
		{
			//while (nTour > (nCluster - (tourNumber * ((nCluster / nUAV))) + 1))
			if (nCluster3 > nUAV)
			{
				float offset = 0.8;
				// selsct boundary tour
				int min_y = 1000000;
				for (int i = 0; i < nTour; i++)
				{
					if (visited[i] == 1);
					else
					{
						for (int j = 0; j < tours[i].size(); j++)
						{
							if (tours[i][j].y < min_y) min_y = tours[i][j].y;
							boundaryTour = i;
						}
					}
				}
				//ns = sizeOfClusters[boundaryTour];
				while ((tours[boundaryTour].size() < offset * (nSensor / nUAV) || tourNumber == nUAV) && (nTour1 > repeat))
				{
					pivot x, minP;
					minP.cost = 1000000;
					for (int i = 0; i < nTour; i++)
					{
						if (i == boundaryTour || visited[i] == 1);
						else
						{
							x = FindMinCycle(tours[boundaryTour], tours[i]);
							if (x.cost < minP.cost)
							{
								jTour = i;
								minP = x;
							}
						}
					}
					//minP.bCross = true;
					tours[boundaryTour] = merge(tours[boundaryTour], tours[jTour], minP);
					//ns += sizeOfClusters[jTour];
					visited[jTour] = 1; visited[boundaryTour] = 1;


					RemoveLoop(&tours[boundaryTour]);
					//while (RemoveCross2(&tours[jTour]));
					nTour1--;

					//cout << "\ncost minP = " << minP.cost << endl;
					final_minP += minP.cost;

					plotSensor << "x, " << "y" << endl;
					//printf("merged tour%d and tour%d\n", jTour, boundaryTour);
					//cout << jTour << endl;

					for (auto i = tours[boundaryTour].begin(); i != tours[boundaryTour].end(); ++i)
					{

						plotSensor << i->x << ", " << i->y << endl;
						//printf("(%d, %d)", i->x, i->y);
					}
					level++;
				}
			}
			visited[boundaryTour] = 1;
			bool f = true;
			while (f)
			{
				f = false;
				while (f = RemoveCross2(&tours[boundaryTour]));
				//printf("\n\n ==> final Cost = %d <==\n", Cost(&tours[jTour]));
				//printf("RemoveCross3----------------\n");

				do
				{
					f = RemoveCross3(&tours[boundaryTour]);
					//printf("\n\n ==> final Cost = %d <==\n", Cost(&tours[jTour]));
				} while (f);
			}


			//printf("\n\n ==> final Cost = %d <==\n", Cost(&tours[boundaryTour]));
			//printf("\n ==> final minP = %d <==\n", final_minP);
			//printf("\n ==> N Sensors  = %d <==\n", tours[boundaryTour].size());
			plotSensor.close();

			//printf("---------Insert loop : begin ------------\n");
			//printf("dist = %d avgAoI = %2.2f\n", Cost(&tours[boundaryTour]), avgAOI(&tours[boundaryTour]));

			while (InsertLoop(&tours[boundaryTour]));

			printf("---------Insert vertex : end ------------\n nSensor = %d\n", tours[boundaryTour].size());
			printf("dist = %d avgAoI = %2.2f\n", Cost(&tours[boundaryTour]), avgAOI(&tours[boundaryTour]));

			while (RemoveCross2(&tours[boundaryTour]));
			while (RemoveCross3(&tours[boundaryTour]));
			printf("after 2/3-opt dist = %d avgAoI = %2.2f\n", Cost(&tours[boundaryTour]), avgAOI(&tours[boundaryTour]));

			ofstream plotTour;
			string tourPlot;
			if (nCluster == 1 && m1 == 2) tourPlot = sensorFileName + ".opt.tour";
			else tourPlot = sensorFileName + " - " + to_string(nCluster) + "cluster - " + method + " - " + to_string(nUAV) + "UAV - " + to_string(tourNumber) + ".tsp.opt.tour";
			plotTour.open(tourPlot);

			ofstream sensorTour;
			string sensorPlot = sensorFileName + " - " + to_string(nCluster) + "cluster - " + to_string(nUAV) + "UAV - " + to_string(tourNumber) + ".csv";
			sensorTour.open(sensorPlot);

			ofstream sensorTsp;
			string sensorT = sensorFileName + " - " + to_string(nCluster) + "cluster - HopCount - " + to_string(nUAV) + "UAV - " + to_string(tourNumber) + ".tsp";
			sensorTsp.open(sensorT);

			ifstream sensorloc;
			string locSensor = sensorFileName;
			sensorloc.open(locSensor);
			for (int i = 0; i < 6; i++) getline(sensorloc, temp);
			int i = 0;
			float x, y;
			int j;
			point s1[SENSOR];
			while (i < nSensor)
			{
				sensorloc >> j >> x >> y;
				s1[i] = point(x, y);
				i++;
			}
			int index_s = 0;
			plotTour << " - " << endl << " - " << endl << " - " << endl << " - " << endl << " - " << endl;
			sensorTsp << " - " << endl << " - " << endl << " - " << endl << " - " << endl << " - " << endl << " - " << endl;
			sensorTour << "x, y" << endl;
			int k = 0;
			bool flag = true;
			int oldIndex = -1;
			for (int i = 0; i < tours[boundaryTour].size(); i++)
			{
				index_s = findSensor(tours[boundaryTour][i], s1);
				if (oldIndex == index_s) continue;
				oldIndex = index_s;
				plotTour << s1[index_s].x << " " << s1[index_s].y << endl;
				sensorTour << s1[index_s].x << ", " << s1[index_s].y << endl;
				flag = true;
				for (int j = 0; j < i; j++) if (index_s == sensorIndexs[j]) flag = false;
				if (flag == true)
				{
					sensorIndexs[k] = index_s;
					k++;
					sensorTsp << s1[index_s].x << " " << s1[index_s].y << endl;
					flag = true;
				}
			}
			sensorTsp << "EOF";
			plotTour << "-1\nEOF";
			sizeOfClusters[tourNumber - 1] = k;

			cout << " ns = " << sizeOfClusters[tourNumber - 1] << endl;
			costs << Cost(&tours[boundaryTour]) << ", " << avgAOI(&tours[boundaryTour]) << ", " << sizeOfClusters[tourNumber - 1] << ", " << tours[boundaryTour].size() << ", ";

			if (nCluster3 == nUAV) boundaryTour++;
			/*if (boundaryTour != 0)
			{
			tours[MAX_TOURS - 1] = tours[0];
			tours[0] = tours[boundaryTour];
			tours[boundaryTour] = tours[MAX_TOURS - 1];
			}
			for (int i = boundaryTour; i < nCluster - 1; i++)
			{
			tours[i] = tours[i + 1];
			}*/
			repeat--;
			nTour1--;
			tourNumber++;

		}
		ofstream clusterSize;
		string clusterSizes = to_string(nSensor) + "Sensors - " + to_string(nCluster) + "Cluster - " + to_string(nUAV) + "UAV - size of clusters.txt";
		clusterSize.open(clusterSizes);
		for (int i = 0; i < nUAV; i++) clusterSize << "Cluster " << (i + 1) << " : " << sizeOfClusters[i] << endl;
		if (m1 != -1) nCluster3 = 10000;
		cout << "##############################################\n##############################################\n\n";
		costs << endl;
	}
	system("pause");
	return 0;
}