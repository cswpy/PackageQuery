#include "pb/util/udeclare.h"

using namespace pb;
int main() {
    //vector<double> data {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    RMatrixXd matrix (3, 3);
    matrix << 1, 4, 7, 2, 5, 8, 3, 6, 9;
    //VectorXd v1 (5);
    //v1 << 1.1, 2.2, 3.3, 4.4, 5.5;
    // Eigen::Map<RMatrixXd> m2 (data.data(), 4, 2);
    //memcpy(&matrix(0, 0), &data[0][0], data.size()*data[0].size()*sizeof(double) );
    //cout << matrix << endl;
    //cout << endl;
    //memcpy(&v1(0), &data[0][0], data.size()*sizeof(double));
    // data.insert(data.end(), v1.begin(), v1.end());
    // cout << data.size() << endl;

    // vector<double> data (v1.begin(), v1.end());

    // // cout << matrix << endl;
    //     for(auto it=data.begin(); it!=data.end(); it++) {
    //     cout << *it << " ";
    // }
    cout << matrix << endl;
    cout << matrix.colwise().sum() / matrix.cols() << endl;
    cout << matrix.cols() << endl;
    return 0;
}