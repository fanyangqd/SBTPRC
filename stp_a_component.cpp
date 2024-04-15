#include"STP_a_component.h"

bool cmp::operator()(vector<double>& a, vector<double>& b) {
	return a[1] > b[1];
};

StpComponent::StpComponent(){
	mRes = 0;
	mCost = 0;
	mWeight = 0;
	mProb = 0.0;
	mFixedCost = 0;
}

StpInstance::StpInstance(){
	mObjVal = 0.0;
	mCpuTime = 0.0;
	mGap = 0.0;
	mState = 0;
	tempTTB = 0.0;
    tempTTB02 = 0.0;
	tempTTB03 = 0.0;
	timeToBest = 0.0;
}

void StpInstance::Initialization() {
	double sumWeight = 0.0;

	for (int i = 1; i <= STP_N; i++){
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(1, 50);
		mComp[i].mCost = dis(gen);

		std::random_device rd02;
		std::mt19937 gen02(rd02());
		std::uniform_int_distribution<> dis02(STP_RES_LB, STP_RES_UB);
		mComp[i].mRes =  dis02(gen02);

		std::random_device rd03;
		std::mt19937 gen03(rd03());
		std::uniform_int_distribution<> dis03(1, 1000);
		mComp[i].mWeight = dis03(gen03);
		sumWeight += double(mComp[i].mWeight);

		mComp[i].mFixedCost = STP_FIXEDCOST;
	}

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(STP_QRANGE_LB, STP_QRANGE_UB);
	double q = dis(gen);
	for (int i = 1; i <= STP_N; i++){
	     mComp[i].mProb = pow(q, double(mComp[i].mWeight) / sumWeight);
	}


	/*// Order components by non-increasing resource requirements 
	for (int i = 1; i <= STP_N; ++i) {
		int max = i;
		for (int j = i + 1; j <= STP_N; ++j) {
			if (this->mComp[j].mRes > this->mComp[max].mRes) {
				max = j;
			}
		}
		StpComponent temp = this->mComp[max];
		this->mComp[max] = this->mComp[i];
		this->mComp[i] = temp;
	}*/

	// Order components by ascending ratio
	for (int i = 1; i <= STP_N; ++i) {
		int min = i;
		for (int j = i + 1; j <= STP_N; ++j) {
			if (this->mComp[j].mCost/(1- this->mComp[j].mProb)< this->mComp[min].mCost/(1-this->mComp[min].mProb)) {
				min = j;
			}
		}
		StpComponent temp = this->mComp[min];
		this->mComp[min] = this->mComp[i];
		this->mComp[i] = temp;
	}


	
	return;
}

void StpInstance::Reset() {
	mObjVal = 0;
	mCpuTime = 0;
	mGap = 0;
	mState = 0;
	for (int i = 1; i <= STP_N; i++)
	{
		mComp[i].mRes = 0;
		mComp[i].mProb = 0.0;
		mComp[i].mCost = 0;
		mComp[i].mWeight = 0;
		mComp[i].mFixedCost = 0;
	}
	return;
}


void StpInstance::Output(int Instance_Number) {

		ofstream outfile;
		ostringstream instanceName;
		instanceName << "STP_" << STP_N << "_Q_" << STP_QRANGE_LB << "_" 
				<< STP_QRANGE_UB << "_Res_" << STP_RES_LB<< "_"<<STP_RES_UB
			<<"_FixedC_"<<STP_FIXEDCOST<<"_"<< Instance_Number << ".txt";
		outfile.open(instanceName.str());

		for (int i = 1; i <= STP_N; ++i) {
			outfile<<fixed << i << " " << mComp[i].mProb << " " << mComp[i].mCost
				<< " " << mComp[i].mRes << " " <<mComp[i].mWeight<<" "<<mComp[i].mFixedCost<<endl;
		}
		outfile.close();
		return;
}



void StpInstance::Input(string Instance_Name) {
	ifstream infile;
	infile.open(Instance_Name);
	this->Reset();
	for (int i = 1; i <= STP_N; ++i) {
		infile >> i >> mComp[i].mProb >> mComp[i].mCost >> mComp[i].mRes >> mComp[i].mWeight >> mComp[i].mFixedCost;
	}
	infile.close();

	return;
}


StpBatch::StpBatch() {

	this->batchCost = 0;
	this->compNumber = 0;
	this->batchRatio = 0.0;
	this->batchProb = 1.0;
	this->batchID = 0;
	this->comps.resize(0);

	return;
}

// Functions:
double max_number(double a, double b) { return (a > b) ? a : b; }

double min_number(double a, double b) {
	if (a < b) { return a; }
	else { return b; }
}
int min_number(int a, int b) {
	if (a < b) { return a; }
	else { return b; }
}

void compute_ObjVal(double& xObj, vector<vector<double>>& xBatchSeq) {
	// Caculate objective value
	sort(xBatchSeq.begin(), xBatchSeq.end(), lessSort);
	double prob = 1.0;
	xObj = 0.0;
	for (int i = 0; i < xBatchSeq.size(); i++) {
		xObj += xBatchSeq[i][3] * prob;
		prob = prob * xBatchSeq[i][2];
	}
	return;
}

bool lessSort(vector<double>& a, vector<double>& b) { return (a[1] < b[1]); }


