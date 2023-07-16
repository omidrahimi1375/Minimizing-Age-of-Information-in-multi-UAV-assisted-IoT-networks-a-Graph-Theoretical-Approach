// AOI.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include "stdlib.h"
#include "string.h"
#include <windows.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <list>
#include <vector>

using namespace std;

#define MAX_X 10000
#define MAX_Y 10000
#define ARRLEN(x) (sizeof(x)/sizeof(x[0]))
#define SENSOR 300
#define EPISOD 20000

int nSensors, nSensors1;
string sensorFileName, sensorLocation;


//#pragma comment(linker, "/STACK:500000000")
//#pragma comment(linker, "/HEAP:500000000")

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

point Action[] = { point(0, 1), point(0, -1), point(1, 0), point(-1, 0) }; // r l d u
enum
{
	Y_PLUS = 0,
	Y_MINUS = 1,
	X_PLUS = 2,
	X_MINUS = 3
};
point Sensors[SENSOR];

struct CheckPoint
{
	int slot_no;
	int sensorID;
	point pos;
	CheckPoint(int slot, int id, point p)
	{
		slot_no = slot;
		sensorID = id;
		pos = p;
	}
	void reset()
	{
		slot_no = -1;
		sensorID = -1;
		pos = point(-1, -1);
	}
};
struct AgeVector
{
	float a[SENSOR];
	int dirUAV;
	int indexUAV;
	AgeVector(int val = 0)
	{
		reset(val);
	}
	void reset(int val)
	{
		for (int i = 0; i < nSensors; i++)
		{
			a[i] = val;
			indexUAV = -1;
			dirUAV = 0;
		}
	}
	float sum(int level = 0)
	{
		float res = 0;
		if (level == 0)
		{
			for (int i = 0; i < nSensors; i++)
				res += a[i];
		}
		else if (level == 1)
		{
			for (int i = 0; i < nSensors; i++)
				res += dirUAV;
		}
		return res;
	}
	int min_index()
	{
		int idx = 0;
		for (int i = 1; i < nSensors; i++)
			if (a[i] < a[idx])
				idx = i;
		return idx;
	}
	AgeVector operator+(AgeVector second)
	{
		AgeVector res;
		for (int i = 0; i < nSensors; i++)
			res.a[i] = a[i] + second.a[i];
		return res;
	}
	AgeVector operator*(float f)
	{
		AgeVector res;
		for (int i = 0; i < nSensors; i++)
			res.a[i] = f * a[i];
		return res;
	}
	char *getstr(char *tmp)
	{
		sprintf(tmp, "%2.2f ", a[0]);
		char t[100];
		for (int i = 1; i < nSensors; i++)
		{
			if (a[i] == 0)
				sprintf(t, "<%d><<%2.2f>> ", i, a[i]);
			else
				sprintf(t, "%2.2f ", a[i]);
			strcat(tmp, t);
		}
		return tmp;
	}
};

struct QTable
{
	int nRow, nCol;
	AgeVector max_res;
public:
	QTable(int n, int m)
	{
		nRow = n;
		nCol = m;
	}
	float find_max(float *temp, int size)
	{
		float max = 0;
		for (int i = 0; i < size; i++)
		{
			if (max < temp[i]) max = temp[i];
		}
		return max;
	}
	void init()
	{
	}
	AgeVector Reward(AgeVector aoi, point pos, int act, int prev_act, bool bPrint = false, vector<CheckPoint> *pVisitedSensors = NULL)
	{
		AgeVector res;
		int t = 0;
		float dFactor = 1;
		for (int i = 0; i < nSensors; i++)
		{
			float d = dFactor*pos.distance(Sensors[i]);
			if (d == 0)
				res.a[i] = 1000;
			else if (pos.IsOnPath(Sensors[i], Action[act]))
			{
				if (aoi.a[i] > d)
					res.a[i] = (-aoi.a[i] + d);
				else
					res.a[i] = (-aoi.a[i]) / (d);
			}
			/*float temp[SENSOR];
			for (int j = 0; j < nSensors; j++)
			{
			point pos2 = Sensors[i];
			float d2 = dFactor*pos2.distance(Sensors[j]);
			if (aoi.a[j] > d2 && j != i)
			temp[j] = (-aoi.a[j] + d2);
			else
			temp[j] = (-aoi.a[j]) / (d2);
			}
			res.a[i] += find_max(temp, nSensors);*/
		}
		if (bPrint)
		{
			printf("UAV : ");
			for (int i = 0; i < ARRLEN(Sensors); i++)
			{
				if (res.a[i] != 0)
					printf(" S%d(%2.2f) ", i, res.a[i]);
			}
			printf("\n");
		}
		AgeVector rew;
		rew.reset(0);
		// check for cycle

		/*if (pVisitedSensors != NULL)
		{
		for (int i = 0; i < pVisitedSensors->size() && i < nSensors - 1; i++)
		{
		// Clearing the aoi of n - 1 previous meeting sensors
		res.a[pVisitedSensors->at(i).sensorID] = 0;
		}
		}*/

		/*
		if (res.a[res.min_index()] == 0)

		printf("Error : no sensor found with valid aoi\n");
		else
		printf("OK : Reward return valid sensor %d \n", res.min_index());
		*/

		rew.a[res.min_index()] = res.a[res.min_index()];
		for (int i = 0; i < nSensors; i++)
			if (res.a[i] > 0)
				rew.a[i] = res.a[i];
		if (act == prev_act)
			rew.dirUAV = -1;
		else if (Action[act] + Action[prev_act] == point(0, 0))
			rew.dirUAV = 1;
		rew.indexUAV = 0;
		return rew;
	}
	bool IsValid(point p)
	{
		if (p.x < 0 || p.x >= nRow || p.y < 0 || p.y >= nCol)
			return false;
		return true;
	}
	void update(AgeVector age, point currentPos, int act)
	{

	}
	int getValMax(AgeVector age, point pos)
	{
		return 0;
	}
	int getMaxIdx(AgeVector age, point pos)
	{
		return 0;
	}

};

QTable Q(6, 6);
//Reward = 0 - 21 0 0 0
//Q(0 38 19 19 0) UAV - 0: Pos = (0, 0), act = <1, 0>  UAV - 1 : Pos = (19, 19), act = <-1, 0> = -21
void convert(int iAction, int nUAV, int res[])
{
	for (int j = 0; j < nUAV; j++)
	{
		res[j] = (iAction >> (2 * j)) & 3;
	}
}
list<int> nSensorsInRect(point, point, bool);
bool initializeSensors()
{
	fstream my_file;
	int nRow = 0, nCol = 0;
	memset(Sensors, 0, SENSOR * sizeof(int));
	my_file.open(sensorFileName, ios::in);
	if (!my_file)
	{
		cout << "\nFile " << sensorFileName << " does not exist.\n";
		return false;
	}
	else
	{
		int num;
		float x, y;
		int i = 0;
		string line;
		for (int i = 0; i < 6; i++) getline(my_file, line);
		while (i < nSensors)
		{
			my_file >> num >> x >> y;
			Sensors[num - 1] = point(x, y);
			i++;
			nRow = __max(nRow, x);
			nCol = __max(nCol, y);
			//printf("nRow = %d, ncol=%d\n", nRow, nCol);
		}
	}
	my_file.close();

	int max_distance = 0;
	int max_x = 0;
	int max_y = 0;
	int x;
	for (int i = 0; i < nSensors; i++)
	{
		if (Sensors[i].x > max_x) max_x = Sensors[i].x;
		if (Sensors[i].y > max_y) max_y = Sensors[i].y;
		for (int j = i + 1; j < nSensors; j++)
		{
			x = Sensors[i].distance(Sensors[j]);
			if (x > max_distance) max_distance = x;
		}
	}
	Q.nRow = nRow + 1;
	Q.nCol = nCol + 1;
	printf("(nRow, nCol) = (%d, %d)\n", Q.nRow, Q.nCol);
	cout << "max distance = " << max_distance << endl;
	cout << "Grid = " << max_x << " * " << max_y << " = " << max_x * max_y << endl;

	return true;
}

void printSensors()
{
	for (int i = 0; i < nSensors; i++)
		cout << "S" << i + 1 << " : " << Sensors[i].x << " - " << Sensors[i].y << endl;
}
void printSensors1()
{
	ofstream plots;
	plots.open(sensorLocation);
	plots << "x, " << "y" << endl;
	for (int i = 0; i < nSensors; i++)
		plots << Sensors[i].x << ", " << Sensors[i].y << endl;
	plots.close();
}

void printAction(int iAction, int nUAV, AgeVector this_rew, AgeVector curr_rew, point curPos[], AgeVector curAOI)
{
	char tmp[100000];
	int action[SENSOR];
	convert(iAction, nUAV, action);
	printf("\nUpdating %f => %f ... \n", curr_rew.sum(), this_rew.sum());
	//printf("prevMaxRew = %s\n\n", curr_rew.getstr(tmp));
	//printf("thisRew = %s\n", this_rew.getstr(tmp));
	for (int i = 0; i < nSensors; i++) int iUAV = this_rew.indexUAV;
}
bool IsBetween(point s, point p1, point p2)
{
	int xLow = __min(p1.x, p2.x);
	int xHigh = __max(p1.x, p2.x);
	int yLow = __min(p1.y, p2.y);
	int yHigh = __max(p1.y, p2.y);
	return (xLow <= s.x && s.x <= xHigh) && (yLow <= s.y && s.y <= yHigh);
}
using namespace std;
list<int> nSensorsInRect(point p1, point p2, bool bPrint = false)
{
	list<int> res;
	int nCount = 0;
	for (int i = 0; i < nSensors; i++)
		if (IsBetween(Sensors[i], p1, p2))
		{
			res.push_back(i);
			nCount++;
			if (bPrint)
				printf("%d => (%d, %d)\n", nCount, Sensors[i].x, Sensors[i].y);
		}
	return res;
}
using namespace std;
#define MAX_UAV 6
struct Path
{
	list<int> paths[MAX_UAV];
	list<int> meetSensors[MAX_UAV];
	int nTotal[MAX_UAV];
	int dist[MAX_UAV];
};
struct PathMetric
{
	int nHop;
	float aoiSum;
};
PathMetric route[MAX_X][MAX_Y];
bool VisitedInPrevMoves(vector<CheckPoint> *pVisitedSensors, int nMove, int sensorID)
{
	for (int i = 1; i < pVisitedSensors->size() && i < nMove; i++)
	{
		if (pVisitedSensors->at(i).sensorID = sensorID)
			return true;
	}
	return false;
}
void FindPath(int nUAV, point curPos[], AgeVector& curAOI, AgeVector& rewAge, int pMetric, Path *res, vector<CheckPoint> *pVisitedSensors = NULL)
{
	bool bLog = false;
	memset(route, 4 * MAX_X*MAX_Y, 0);
	point dst[SENSOR] = { (-1, -1) };
	dst[0] = Sensors[rewAge.min_index()];
	if (rewAge.min_index() < 0 || rewAge.min_index() >= nSensors)
	{
		printf("Error : destination is not found.\n");
		exit(1);
	}

	for (int iUAV = 0; iUAV < nUAV; iUAV++)
	{
		list<int> sn = nSensorsInRect(curPos[iUAV], dst[iUAV]);
		bool bPrint = (sn.size() > 2);
		int n1 = sn.size();
		if (bLog && bPrint)
		{
			printf("no.%d : %d Sensors in Rect of UAV-%d\n", pVisitedSensors->size(), n1, iUAV);
			while (sn.size())
			{
				int nS = sn.front();
				sn.pop_front();
				printf("S[%d] : %2.2f \n", nS, curAOI.a[nS]);
			}
		}
		//printf("FindPath Pu%d=(%d, %d) Ps=(%d, %d) nSensorsInRect=%d\n", 
		//	iUAV, curPos[iUAV].x, curPos[iUAV].y, dst[iUAV].x, dst[iUAV].y, n1);
		res->paths[iUAV].clear();
		int dx = dst[iUAV].x - curPos[iUAV].x;
		int dy = dst[iUAV].y - curPos[iUAV].y;
		int vx = (dx > 0) ? X_PLUS : X_MINUS;
		int vy = (dy > 0) ? Y_PLUS : Y_MINUS;
		int sx = (dx > 0) ? 1 : -1;
		int sy = (dy > 0) ? 1 : -1;
		if (abs(dx) >= MAX_X || abs(dy) >= MAX_Y)
		{
			printf("Error , MAX_X is small ...\n Exiting Program...\n");
			exit(1);
		}
		res->nTotal[iUAV] = nSensorsInRect(curPos[iUAV], dst[iUAV]).size();
		res->dist[iUAV] = abs(dx) + abs(dy);

		for (int i = 0; i <= abs(dx); i++)
		{
			for (int j = 0; j <= abs(dy); j++)
			{
				point p = dst[iUAV];
				p.x += -i * sx;
				p.y += -j * sy;
				list<int> sn = nSensorsInRect(p, p);
				int isSensor = sn.size();
				//if (sn.size() && VisitedInPrevMoves(pVisitedSensors, 10, sn.front()))
				//	isSensor = 0;
				// this sensor must be not in the path of prev UAVs
				float senAOI = (isSensor) ? curAOI.a[sn.front()] : 0;

				if (i == 0 && j == 0)
				{
					route[i][j].nHop = 1;
					route[i][j].aoiSum = senAOI;
				}
				else if (i == 0)
				{
					route[i][j].nHop = isSensor + route[i][j - 1].nHop;
					route[i][j].aoiSum = senAOI + route[i][j - 1].aoiSum;
				}
				else if (j == 0)
				{
					route[i][j].nHop = isSensor + route[i - 1][j].nHop;
					route[i][j].aoiSum = senAOI + route[i - 1][j].aoiSum;
				}
				else
				{
					route[i][j].nHop = isSensor + __max(route[i - 1][j].nHop, route[i][j - 1].nHop);
					route[i][j].aoiSum = senAOI + __max(route[i - 1][j].aoiSum, route[i][j - 1].aoiSum);
				}
			}

		}
		int i = abs(dx);
		int j = abs(dy);
		//printf(" nAvailable=%d\n", route[i][j]);
		res->meetSensors[iUAV].clear();
		while (i > 0 || j > 0)
		{
			point p = dst[iUAV];
			p.x += -i * sx;
			p.y += -j * sy;
			list<int> sn = nSensorsInRect(p, p);
			if (sn.size() > 0)
				res->meetSensors[iUAV].push_back(sn.front());
			if (i > 0 && j > 0)
			{
				float x1 = (pMetric == 1) ? route[i - 1][j].nHop : route[i - 1][j].aoiSum;
				float x2 = (pMetric == 1) ? route[i][j - 1].nHop : route[i][j - 1].aoiSum;
				if (x1 >= x2) // >= bishtarin node ro rad kone
				{
					res->paths[iUAV].push_back(vx);
					i--;
				}
				else
				{
					res->paths[iUAV].push_back(vy);
					j--;
				}
			}
			else if (i > 0)
			{
				res->paths[iUAV].push_back(vx);
				i--;
			}
			else
			{
				res->paths[iUAV].push_back(vy);
				j--;
			}
		}
		if (bLog && res->meetSensors[iUAV].size() > 1)
		{
			printf("Sensors are passing through by UAV-%d => {", iUAV);
			for (int x : res->meetSensors[iUAV])
				if (curAOI.a[x] != 0)
					printf(" S[%d] : %2.2f ", x, curAOI.a[x]);
			printf("}\n");
		}
	}
}

struct Cycle
{
	bool found;
	int weight;
	int nHop;
	int sID;
	int slot_begin, slot_end;
	float AvgAOI;
	int nVisited;
	Cycle()
	{
		weight = nHop = 0;
		nVisited = 0;
	}
	void reset()
	{
		found = 0;
		weight = 0;
		nHop = 0;
		sID = 0;
		slot_begin = 0;
		slot_end = 0;
		AvgAOI = 0;
		nVisited = 0;
	}
};
int getVisitedCount(int Visited[])
{
	int n = 0;
	for (int i = 0; i < nSensors; i++)
		if (Visited[i] == 1)
			n++;
	return n;
}
Cycle FindCycle(vector<CheckPoint> *uP, vector<AgeVector> *sAOI)
{
	CheckPoint first = uP->front();
	int Visited[SENSOR];
	memset(Visited, 0, sizeof(int)*nSensors);
	Cycle res;
	int i;
	for (i = 1; i < uP->size(); i++)
	{
		Visited[uP->at(i).sensorID] = 1;
		res.weight += uP->at(i).pos.distance(uP->at(i - 1).pos);
		if (getVisitedCount(Visited) == nSensors &&
			first.sensorID == uP->at(i).sensorID)
			break;
	}
	res.sID = first.sensorID;
	res.slot_end = first.slot_no;
	res.nHop = i;
	if (i < uP->size())
	{
		res.slot_begin = uP->at(i).slot_no;
	}
	else
	{
		res.slot_begin = uP->at(i - 1).slot_no;
	}
	res.nVisited = getVisitedCount(Visited);
	res.found = (res.nVisited == nSensors);
	if (res.found)
	{
		res.AvgAOI = 0.0;
		for (int i = res.slot_begin; i <= res.slot_end; i++)
			res.AvgAOI += sAOI->at(i).sum() / nSensors;
		res.AvgAOI = res.AvgAOI / (res.slot_end - res.slot_begin);
	}
	return res;
}

float euclidean_dist(point p1, point p2)
{
	return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2) * 1.0);
}

float manhatan_dist(point p1, point p2)
{
	return (fabs(p2.x - p1.x) + fabs(p2.y - p1.y));
}
point findUavPos(int s)
{
	point pos = { Sensors[0].x, Sensors[0].y };
	switch (s)
	{
	case 1:
		for (int i = 1; i < nSensors; i++)
		{
			if (Sensors[i].x <= pos.x && Sensors[i].y >= pos.y)
			{
				pos.x = Sensors[i].x;
				pos.y = Sensors[i].y;
			}
		}
		break;
	case 2:
		for (int i = 1; i < nSensors; i++)
		{
			if (Sensors[i].x <= pos.x && Sensors[i].y <= pos.y)
			{
				pos.x = Sensors[i].x;
				pos.y = Sensors[i].y;
			}
		}
		break;
	case 3:
		for (int i = 1; i < nSensors; i++)
		{
			if (Sensors[i].x >= pos.x && Sensors[i].y >= pos.y)
			{
				pos.x = Sensors[i].x;
				pos.y = Sensors[i].y;
			}
		}
		break;
	case 4:
		for (int i = 1; i < nSensors; i++)
		{
			if (Sensors[i].x >= pos.x && Sensors[i].y <= pos.y)
			{
				pos.x = Sensors[i].x;
				pos.y = Sensors[i].y;
			}
		}
		break;

	}
	return pos;
}
int findSensor(point a)
{
	int finded = -1;
	for (int i = 0; i < nSensors; i++)
	{
		if (a.x == Sensors[i].x && a.y == Sensors[i].y)
		{
			finded = i;
		}
	}
	return finded;
}
int episod = 0;
int main()
{
	int nClusters = 2;

	int clusterFlag = 1;

	while (nClusters < 21 && clusterFlag == 1)
	{
		LARGE_INTEGER tBegin, tEnd, tFreq;
		float final_tDiff = 0;
		float final_cost = 0;


		AgeVector curAOI;
		Path res;
		int mode = 1;
		int nClusters2 = 1;
		string fileName, temp;
		int pathMetric = 2;
		ifstream config;
		config.open("config.txt");
		config >> temp >> temp >> temp;
		sensorFileName = temp + ".tsp";
		sensorLocation = temp + ".csv";
		config >> temp >> temp >> temp;
		nSensors = stoi(temp);
		nSensors1 = nSensors;
		config >> temp >> temp >> temp;
		if (temp == "no") clusterFlag = 0; else clusterFlag = 1;
		config >> temp >> temp >> temp;
		if (temp == "ourMethod") mode = 1; else mode = 2;
		config >> temp >> temp >> temp;
		if (temp == "averageAOI") pathMetric = 2; else if (temp == "hopCount") pathMetric = 1; else pathMetric = -1;
		config >> temp >> temp >> temp;
		int initAOI = stoi(temp);
		curAOI.reset(initAOI);
		config >> temp >> temp >> temp;
		int uavPos;
		if (temp == "LU") uavPos = 1;
		else if (temp == "LD") uavPos = 2;
		else if (temp == "RU") uavPos = 3;
		else uavPos = 4;
		point curPos[1];
		if (clusterFlag == 1)
		{
			config >> temp >> temp >> temp;
			nClusters2 = stoi(temp);
		}
		else nClusters = 1;
		//point curPos[] = { Sensors[stoi(temp)] };
		//point curPos[] = { point(650, 1000) };


		/*cout << "Number of Sensors : \n";
		cin >> nSensors;
		nSensors1 = nSensors;
		cout << "Clustering = 1, else 0 : \n";
		cin >> clusterFlag;
		cout << "Choose mode (1 or 2) : \nOur Method = 1, Optimal Method = 2\n";
		cin >> mode;
		cout << "Select Metric for Path (1 : HopCount, 2: AOISum)\n";
		cin >> pathMetric;
		if (mode == 2)
		{
		cout << "Choose method (1 or 2) : \nOur Method = 1, Optimal Method = 2\n";
		cin >> opt_method;

		}
		if (clusterFlag == 1)
		{
		cout << "Number of clusters : \n";
		cin >> nClusters;
		}
		else {
		switch (nSensors)
		{
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
		}*/


		// UAV Position : {(39, 26), (62, 26), (60, 39), (31, 38)}
		//int aoi_v0[] = { 107,25,214,73,273,254,72,62,208,136,112,194,21,43,293,216,250,168,34,285,252,234,134,153,161,260,276,222,51,60,132,228,242,83,57,165,297,129,213,188,92,62,86,211,129,96,144,42,35,174,205,89,50,9,159,126,10,92,39,181,211,334,223,71,122,102,102,243,155,279,175,194,236,324,112,109 };



		int repeat = 1;
		cout << endl << endl;
		ifstream clusterSizes;
		string t1;
		if (clusterFlag == 1)
		{
			string s = to_string(nSensors) + "Sensors - " + to_string(nClusters) + "Cluster - size of clusters.txt";
			clusterSizes.open(s);
		}
		while (repeat <= nClusters)
		{
			vector<CheckPoint> uPath;
			vector<CheckPoint> uTargetedSensors;

			if (clusterFlag == 1)
			{
				fileName = to_string(nSensors1) + "Sensors - " + to_string(nClusters) + "Cluster - " + to_string(repeat);
				sensorFileName = fileName + ".tsp";
				sensorLocation = fileName + ".csv";
				cout << "Cluster " << repeat << " : \n";
				clusterSizes >> t1 >> t1 >> t1;
				clusterSizes >> nSensors;
				cout << "Size of cluster = " << nSensors << endl;
			}
			if (initializeSensors() == false)
			{
				cout << " false sensors";
				system("pause");
				return 1;
			}
			point pos = findUavPos(uavPos);
			curPos[0].x = 8; curPos[0].y = 6;
			cout << "curPos = " << curPos[0].x << ", " << curPos[0].y << endl;

			//curPos[0] = point(33, 33);
			printSensors1();

			//printSensors();
			//printf("Press Enter ...\n");
			//getchar();

			char tmp[100000];

			// point(1, 1), point(2, 1), point(2, 2), point(2, 3)


			vector<AgeVector> sAOI;

			int nUAV = ARRLEN(curPos);
			/*for (int i = 0; i < ARRLEN(aoi_v0); i++)
			curAOI.a[i] = aoi_v0[i];
			*/

			//curAOI.a[0] = 8;
			//curAOI.a[1] = 11;
			//curAOI.a[2] = 14;
			//curAOI.a[3] = 6;
			//curAOI.a[4] = 9;


			AgeVector one(1);
			Q.init();
			//Q.print(curAOI);

			int nState = (1 << 2 * nUAV);
			int iBestAction = 0;

			if (mode == 1)
			{
				ofstream plots;
				string plot;
				ofstream tours;
				string tour;
				if (pathMetric == 1)
				{
					plot = "plots - ";
					plot += sensorFileName;
					plot += " - HopCount.csv";
				}
				else
				{
					plot = sensorFileName;
					plot += " - AOISum.csv";
				}
				plots.open(plot);

				if (pathMetric == 1)
				{
					tour = sensorFileName;
					tour += " - HopCount.opt.tour";
				}
				else
				{
					tour = sensorFileName;
					tour += " - AOISum.opt.tour";
				}
				tours.open(tour);
				if (nSensors == 1)
				{
					tours << "test\ntest\ntest\ntest\ntest\n";
					tours << Sensors[0].x << " " << Sensors[0].y << endl;

					tours << -1 << endl << "EOF";
					tours.close();
				}
				else
				{
					ofstream route;
					string r = "route - ";
					r += sensorFileName;
					r += ".csv";
					route.open(r);

					route << "route" << endl;

					for (int i = 0; i < nUAV; i++) plots << "x-uav" << i + 1 << ", " << "y-uav" << i + 1 << ", ";
					plots << "Average AOI - " << nUAV << " UAVs - " << nSensors << " Sensors, ";
					for (int i = 0; i < nSensors; i++) plots << "AOI S" << i + 1 << ", ";
					plots << endl;
					for (int j = 0; j < nUAV; j++)
					{
						plots << curPos[j].x << ", " << curPos[j].y << ", ";
					}
					plots << (curAOI.sum() / nSensors) << ", ";
					for (int j = 0; j < nSensors; j++)
					{
						plots << curAOI.a[j] << ", ";
					}
					plots << endl;
					int nMeet = 0;
					bool bPathExist = false;
					bool bMeet = false;

					QueryPerformanceCounter(&tBegin);
					QueryPerformanceFrequency(&tFreq);

					point uavPos[EPISOD];

					Cycle minCycle;
					minCycle.reset();
					minCycle.weight = 1000000;
					curAOI.reset(initAOI);

					ofstream dist;
					dist.open("eil76-dist.csv");
					dist << "episod\n dist\n";

					for (episod = 0; episod < EPISOD; episod++)
					{
						uavPos[episod] = curPos[0];
						int action[SENSOR];
						int prev_action[SENSOR] = { 0 };
						AgeVector max_rew;
						for (int i = 0; !bPathExist && i < nState; i++)
						{
							convert(i, nUAV, action);
							bool isValid = true;
							for (int k = 0; k < nUAV; k++)
							{
								point next1 = curPos[k] + Action[action[k]];
								if (Q.IsValid(next1) == false)
								{
									isValid = false;
									break;
								}
							}
							if (isValid == false) continue;
							AgeVector rewAge = Q.Reward(curAOI, curPos[0], action[0], prev_action[0], false, &uPath);
							if (rewAge.sum() < max_rew.sum())
							{
								iBestAction = i;
								max_rew = rewAge;
							}
							else if (rewAge.sum() == max_rew.sum())
							{
								if (rewAge.sum(1) < max_rew.sum(1))
								{
									iBestAction = i;
									max_rew = rewAge;
								}
								else if (rewAge.sum(1) == max_rew.sum(1) && rand() % 2 == 0)
								{
									iBestAction = i;
									max_rew = rewAge;
								}
							}
						}
						convert(iBestAction, nUAV, action);
						convert(iBestAction, nUAV, prev_action);

						AgeVector rewAge = Q.Reward(curAOI, curPos[0], action[0], prev_action[0], false, &uPath);
						if (bPathExist == false)
						{
							nMeet = 0;
							FindPath(nUAV, curPos, curAOI, rewAge, pathMetric, &res, &uPath);
							uTargetedSensors.insert(uTargetedSensors.begin(), CheckPoint(episod, rewAge.min_index(), curPos[0]));
							if (res.paths[0].size() == 0)
							{
								printf("Error in path finding ... \n");
								rewAge = Q.Reward(curAOI, curPos[0], action[0], prev_action[0], false, &uPath);
								int g = 1;
							}
							bPathExist = true;
						}
						for (int j = 0; j < nUAV; j++)
						{
							action[j] = res.paths[j].front();
							res.paths[j].pop_front();
							if (res.paths[j].size() == 0)
								bPathExist = false;
							curPos[j] = curPos[j] + Action[action[j]];
						}
						curAOI = curAOI + one;
						bMeet = false;
						for (int i = 0; i < nSensors; i++)
						{
							for (int j = 0; j < nUAV; j++)
							{
								if (curPos[j] == Sensors[i])
								{
									curAOI.a[i] = 0;
									nMeet++;
									uPath.insert(uPath.begin(), CheckPoint(episod, i, curPos[j]));
									bMeet = true;
									char tmp[4000];
								}
							}
						}
						sAOI.push_back(curAOI);
						if (bMeet)
						{
							Cycle cy = FindCycle(&uPath, &sAOI);
							if (cy.found && (cy.weight <= minCycle.weight)/* || (cy.AvgAOI  < minCycle.AvgAOI )*/)
							{
								minCycle = cy;
								printf("Finding Cycle sID = %d weight = %d nHop = %d slot:%d->%d avgAOI = %2.2f\n",
									cy.sID, cy.weight, cy.nHop, cy.slot_begin, cy.slot_end, cy.AvgAOI);
								dist << cy.slot_end << ", " << cy.weight << endl;
							}
						}
						for (int j = 0; j < nUAV; j++)
						{
							plots << curPos[j].x << ", " << curPos[j].y << ", ";
						}
						plots << (curAOI.sum() / nSensors) << ", ";
						for (int j = 0; j < nSensors; j++)
						{
							plots << curAOI.a[j] << ", ";
						}
						plots << endl;
					}
					plots.close();

					Cycle cy = minCycle;
					printf("Min Cycle sID = %d weight = %d nHop = %d slot:%d->%d avgAOI = %2.2f\n",
						cy.sID, cy.weight, cy.nHop, cy.slot_begin, cy.slot_end, cy.AvgAOI);

					string uavRoute;
					if (pathMetric == 1)
					{
						uavRoute = "route - ";
						uavRoute += sensorFileName;
						uavRoute += " - HopCount - ";
						uavRoute += to_string(nUAV);
						uavRoute += "UAV.opt.csv";
					}
					else if (pathMetric == 2)
					{
						uavRoute = "route - ";
						uavRoute += sensorFileName;
						uavRoute += " - AOISum - ";
						uavRoute += to_string(nUAV);
						uavRoute += "UAV.opt.csv";
					}
					ofstream ur;
					ur.open(uavRoute);
					ur << "x, " << "y" << endl;
					tours << "test\ntest\ntest\ntest\ntest\n";
					for (int i = minCycle.slot_begin; i <= minCycle.slot_end; i++)
					{
						for (int j = 0; j < nSensors; j++)
						{
							if (uavPos[i] == Sensors[j])
							{
								tours << Sensors[j].x << " " << Sensors[j].y << endl;
								ur << Sensors[j].x << ", " << Sensors[j].y << endl;
							}
						}
					}
					tours << -1 << endl << "EOF";
					tours.close();
					route.close();
				}
			}

			else if (mode == 2)
			{
				ofstream plots;
				string plot;
				if (pathMetric == 1)
				{
					plot = "plots - ";
					plot += sensorFileName;
					plot += " - HopCount - ";
					plot += to_string(nUAV);
					plot += "UAV.opt.csv";
				}
				else if (pathMetric == 2)
				{
					plot = "plots - ";
					plot += sensorFileName;
					plot += " - AOISum - ";
					plot += to_string(nUAV);
					plot += "UAV.opt.csv";
				}
				else
				{
					plot = "plots - ";
					plot += sensorFileName;
					plot += " - ";
					plot += to_string(nUAV);
					plot += "UAV.opt.csv";
				}
				plots.open(plot);

				for (int i = 0; i < nUAV; i++) plots << "x-uav" << i + 1 << ", " << "y-uav" << i + 1 << ", ";
				plots << "Average AOI - " << nUAV << " UAVs - " << nSensors << " Sensors, ";
				for (int i = 0; i < nSensors; i++) plots << "AOI S" << i + 1 << ", ";
				plots << endl;

				ifstream opt_tour;
				string optimal_tour;
				if (pathMetric == 1)
				{
					optimal_tour = sensorFileName;
					optimal_tour += " - HopCount.opt.tour";
				}
				else if (pathMetric == 2)
				{
					optimal_tour = sensorFileName;
					optimal_tour += " - AOISum.opt.tour";
				}
				else
				{
					optimal_tour = sensorFileName;
					optimal_tour += ".opt.tour";
				}
				opt_tour.open(optimal_tour);
				int index = -1;

				int i = 0;
				string line;
				vector <point> UAV_pos_opt;
				vector <point> Sensors_opt;
				//vector<point> UAV_pos_opt;
				//vector<int> Sensors_opt;

				int x1, y1;
				for (int i = 0; i < 5; i++) getline(opt_tour, line);
				/*
				ofstream opt1;
				opt1.open("a280-1.txt");
				opt1 << " - \n - \n - \n - \n - \n";
				int x = 0;
				while (x != -1)
				{
				opt_tour >> x;
				if(x != -1)
				opt1 << Sensors[x - 1].x << " " << Sensors[x - 1].y << endl;
				}
				*/
				int x0, y0;
				opt_tour >> x1 >> y1;
				x0 = x1; y0 = y1;
				curPos[0].x = x1; curPos[0].y = y1;
				QueryPerformanceCounter(&tBegin);
				QueryPerformanceFrequency(&tFreq);
				bool end = false;
				while (nSensors != 1)
				{
					int x = x1 - curPos[0].x;
					int y = y1 - curPos[0].y;

					if (x == 0 && y == 0)
					{
						if (UAV_pos_opt.size() != 0 && x1 == x0 && y1 == y0 && end == true) break;
						Sensors_opt.push_back(point(x1, y1));
						opt_tour >> x1 >> y1;

					}
					else
					{
						UAV_pos_opt.push_back(curPos[0]);
						if (x > 0) curPos[0].x++;
						else if (x < 0) curPos[0].x--;
						else if (y > 0) curPos[0].y++;
						else if (y < 0) curPos[0].y--;

					}
					if (x1 == -1)
					{
						x1 = x0; y1 = y0, end = true;
					}
				}


				float cost = 0;
				int v = Sensors_opt.size();
				cout << "Size of cluster = " << v << endl;
				int u = UAV_pos_opt.size();
				for (int i = 0; i < v; i++)
				{
					cost += manhatan_dist(Sensors_opt[i], Sensors_opt[(i + 1) % v]);
					//cost += Sensors[Sensors_opt[i]].distance(Sensors[Sensors_opt[(i + 1) % v]]);
				}
				cout << "cost = " << cost << endl;
				final_cost += cost;

				int s[10];
				for (int i = 0; i < nUAV; i++) s[i] = (u / nUAV) * i;

				int k = 0;

				//for (int i = 0; i < nUAV; i++) k[i] = (v / nUAV) * i;

				int RNDMAX = 10;
				char str[10000];
				for (int rnd = 0; rnd < RNDMAX; rnd++)
				{
					//curAOI.getstr(str);
					//printf("curAOI = %s\n", str);
					float sum = 0;
					for (int i = 0; i < u; i++)
					{
						curAOI = curAOI + one;
						curPos[0] = UAV_pos_opt[i];

						int finded = findSensor(curPos[0]);
						if (finded != -1)
						{
							curAOI.a[finded] = 0;
						}
						sum += curAOI.sum() / nSensors;
						plots << curPos[0].x << ", " << curPos[0].y << ", ";
						plots << (curAOI.sum() / nSensors) << ", ";
						for (int j = 0; j < nSensors; j++)
						{
							plots << curAOI.a[j] << ", ";
						}
						plots << endl;
					}
					printf("rnd = %d aoi = %4.4f aoi_avg = %f\n", rnd, curAOI.sum() / nSensors, sum / u);
					//curAOI.getstr(str);
					//printf("curAOI = %s\n", str);
					//getchar();
					//for (int j = 0; j < nUAV; j++)
					//{

					//}

				}
				//for (int i = 0; i < u; i++) cout << i << " : " << UAV_pos_opt[i].x << ", " << UAV_pos_opt[i].y << endl;
				plots.close();

				string route;
				if (pathMetric == 1)
				{
					route = "route - ";
					route += sensorFileName;
					route += " - HopCount - ";
					route += to_string(nUAV);
					route += "UAV.opt.csv";
				}
				else if (pathMetric == 2)
				{
					route = "route - ";
					route += sensorFileName;
					route += " - AOISum - ";
					route += to_string(nUAV);
					route += "UAV.opt.csv";
				}
				else
				{
					route = "route - ";
					route += sensorFileName;
					route += " - ";
					route += to_string(nUAV);
					route += "UAV.opt.csv";
				}
				ofstream r;
				r.open(route);
				r << "x, " << "y" << endl;
				for (int i = 0; i < v; i++)
					r << Sensors_opt[i].x << ", " << Sensors_opt[i].y << endl;
				r.close();
			}

			QueryPerformanceCounter(&tEnd);
			float tDiff = (float)(tEnd.QuadPart - tBegin.QuadPart) / tFreq.QuadPart;
			printf(" ==> Time : %2.4f <==\n", tDiff);
			final_tDiff += tDiff;

			repeat++;
			cout << endl;
		}
		printf("\n ==> Final Time : %2.4f <==\n", final_tDiff);
		printf("\n ==> Final Cost : %2.4f <==\n\n", final_cost);
		/*printf("\n\nVisited Sensors\n\n");
		int nCount = 0;
		for (int i = uPath.size() - 1; i >= 0; i--)
		{
		// Clearing the aoi of n - 1 previous meeting sensors
		printf("%3d ", uPath.at(i).sensorID);
		if (nCount++ % nSensors == nSensors - 1)
		printf("\n");
		}
		nCount = 0;
		printf("\n\nTageted Sensors\n\n");
		for (int i = uTargetedSensors.size() - 1; i >= 0; i--)
		{
		// Clearing the aoi of n - 1 previous meeting sensors
		printf("%3d ", uTargetedSensors.at(i).sensorID);
		if (nCount++ % nSensors == nSensors - 1)
		printf("\n");
		}*/
		nClusters++;
	}
	system("pause");
	return 0;

}