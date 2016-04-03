#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <QMainWindow>
#include <string>
#include <QtGnuplot/QtGnuplotWidget.h>
#include <QtGnuplot/QtGnuplotInstance.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <lapacke.h>
#include <iostream>

#include "cpp/initfuncs.h"
#include "cpp/polyfuncs.h"
#include "cpp/sylvester.h"
#include "cpp/polys.h"
#include "cpp/sy.h"
#include "cpp/compan.h"
#include "cpp/generalized.h"
#include "cpp/interpolation.h"

#define THRESHOLD 0.00001

using namespace Eigen;
using namespace std;

namespace Ui {
class Visualization;
}

class Visualization : public QMainWindow
{
    Q_OBJECT

public:
    explicit Visualization(QWidget *parent = 0);
    ~Visualization();
    void printall();

    void solve();
    void file_read();

    void generate_interpolation();
    bool isready() { return flag; }

    string get_func1(void){ return func1; }
    string get_func2(void){ return func2; }

private:
    QtGnuplotWidget *widget;
    QtGnuplotInstance *instance;

    ofstream output;

    int flagk;
    int k1, k2;
    int counter;
    int roots_count;
    int deg1, deg2;
    double max_x,max_y;
    bool flag;
    bool varchange;
    char main, hidden;
    string func1,func2;
    BivPoly *f1;
    BivPoly *f2;
    MatrixXd pts1, pts2;
    QString roots;
    QString xrange, yrange;
    QString sol_xrange, sol_yrange;
    QStringList points1, points2;

protected:
    bool eventFilter(QObject *obj, QEvent *event);

private slots:
    void on_actionExit_triggered();

    void on_pushButton_clicked();

    void on_pushButton_3_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_7_clicked();

    void on_pushButton_8_clicked();

    void on_pushButton_9_clicked();

    void on_pushButton_2_clicked();

    void on_checkBox_toggled(bool checked);

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

private:
    Ui::Visualization *ui;
};

/******* Help Functions *******/
MatrixXd qlist_to_matrix(QStringList lst);
QString matrix_to_qstring(MatrixXd mtrx);

#endif // VISUALIZATION_H
