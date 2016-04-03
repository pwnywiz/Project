#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>
#include <QString>
#include <stdlib.h>
#include <math.h>
#include "visualization.h"
#include "ui_visualization.h"

using namespace std;

//Constructor for Visualization
Visualization::Visualization(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Visualization)
{
    ui->setupUi(this);
    output.open("roots");
    flag = 0;
    flagk = 0;
    deg1 = -1;
    deg2 = -1;
    counter = 0;
    roots_count = 0;
    varchange = false;
    max_x = 0;
    max_y = 0;
    widget = new QtGnuplotWidget();
    widget->installEventFilter(this);
    widget->setStatusLabelActive(true);
    instance = new QtGnuplotInstance();
    instance->setWidget(widget);
}

/*****************************************************************************
******************************************************************************/
//Destructor for Visualization
Visualization::~Visualization()
{
    delete ui;
    delete instance;
    delete widget;
    delete f1;
    delete f2;
    output.close();
}

/*****************************************************************************
******************************************************************************/
//menu bar,exit option
void Visualization::on_actionExit_triggered()
{
    qApp->quit();
}

/*****************************************************************************
******************************************************************************/
//Setting points in plot widget
//widget auto-closes when k points were given
//and prints them in the correct listwidget
bool Visualization::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::MouseButtonPress)
    {
        if (obj == this->widget) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent *>(event);
            if (mouseEvent->button() == Qt::LeftButton) {
                counter++;
                if (flagk == 1) {
                    points1 << this->widget->getStatusLabel()->text();
                    if (counter == k1) {
                        counter = 0;
                        this->widget->hide();
                        ui->listWidget_1->clear();
                        ui->listWidget_1->addItems( points1 );
                    }
                }
                if (flagk == 2) {
                    points2 << this->widget->getStatusLabel()->text();
                    if (counter == k2) {
                        counter = 0;
                        this->widget->hide();
                        ui->listWidget_2->clear();
                        ui->listWidget_2->addItems( points2 );
                    }
                }
            }
        }
    }
    return QObject::eventFilter(obj, event);
}

/*****************************************************************************
******************************************************************************/
//Reading equations from given file
//Polynomial's analysis algorithm
void Visualization::file_read() {
    int mncounter1 = 1, mncounter2 = 1;

    int maxexpv;

    bool checkdeg1(true);
    bool checkdeg2(true);


    int x1,y1,x2,y2;

    mncounter1=countmononyms(func1);
    mncounter2=countmononyms(func2);

    //buffer arrays for classes creation
    int* expx1=new int[mncounter1];
    int* expy1=new int[mncounter1];
    double* coeff1=new double[mncounter1];

    int* expx2=new int[mncounter2];
    int* expy2=new int[mncounter2];
    double* coeff2=new double[mncounter2];

    checkdeg1 = polanalysis(func1,deg1,mncounter1,expx1,expy1,coeff1);
    if(!checkdeg1){
        ui->textBrowser->append("The given degree of 1st function was wrong and fixed to: "+QString::number(deg1)+"\n");
    }

    x1=findmax(expx1,mncounter1);
    y1=findmax(expy1,mncounter1);

    checkdeg2 = polanalysis(func2,deg2,mncounter2,expx2,expy2,coeff2);
    if(!checkdeg2){
        ui->textBrowser->append("The given degree of 2nd function was wrong and fixed to: "+QString::number(deg2)+"\n");
    }

    x2=findmax(expx2,mncounter2);
    y2=findmax(expy2,mncounter2);

    maxexpv=x1;
    main='x';
    hidden='y';

    if(x2>maxexpv){
        maxexpv=x2;
    }
    if(y1>maxexpv){
        maxexpv=y1;
        main='y'; hidden='x';
    }
    if(y2>maxexpv){
        maxexpv=y2;
        main='y'; hidden='x';
    }

    if(main=='x'){
        f1 = new BivPoly(expx1,expy1,coeff1,x1,y1,mncounter1);
        f2 = new BivPoly(expx2,expy2,coeff2,x2,y2,mncounter2);
    }
    else{
        f1 = new BivPoly(expy1,expx1,coeff1,y1,x1,mncounter1);
        f2 = new BivPoly(expy2,expx2,coeff2,y2,x2,mncounter2);
    }

    //delete temp arrays for class creation,we dont need them anymore
    delete[] expx1;
    delete[] expy1;
    delete[] coeff1;
    delete[] expx2;
    delete[] expy2;
    delete[] coeff2;
}

/*****************************************************************************
******************************************************************************/
//Interpolation computation function
void Visualization::generate_interpolation() {

    pts1 = qlist_to_matrix(points1);
    pts2 = qlist_to_matrix(points2);

    //Uncomment the following in case of hard coded example
    //the two lines before must be commented !

//    pts1.resize(5,2);
//    pts2.resize(5,2);


    /*EXAMPLE 1*/
//    pts1 << -1,0,
//          4,-1,
//          -1,-5,
//          -4,2,
//          4,-4;

//    pts2 << -1,-5,
//          -3,3,
//          -4,-4,
//          -3,-1,
//          -2,-2;


    /*EXAMPLE 2*/
//    pts1 << -1,-5,
//          -1,0,
//          -3,2,
//           3,3,
//           -1,-3;

//    pts2 << -5,3,
//          2,4,
//          5,1,
//          -3,1,
//          -3,2;

    /*EXAMPLE 3*/
//    pts1 << -5,1,
//          -5,0,
//          2,-2,
//          -5,-2,
//           4,4;

//    pts2 << -3,-4,
//          0,4,
//          -5,0,
//          1,-5,
//          -4,3;

    Interpolation inter1(pts1,deg1);
    Interpolation inter2(pts2,deg2);
    QString in1 = matrix_to_qstring(inter1.get_inter());
    ui->textBrowser->append("1st Interpolation =");
    ui->textBrowser->append(in1);
    QString in2 = matrix_to_qstring(inter2.get_inter());
    ui->textBrowser->append("2nd Interpolation =");
    ui->textBrowser->append(in2);

    if (!inter1.solver() || !inter2.solver()) {
      ui->textBrowser->append("Stopping the solution of functions...\n");
      exit(1);
    }

    //prints out interpolation-related information tou App.output
    ui->textBrowser->append("1st Rank = "+QString::number(inter1.get_rank())+"\n");
    QString ker1 = matrix_to_qstring(inter1.get_kernel());
    ui->textBrowser->append("1st Kernel =");
    ui->textBrowser->append(ker1);
    ui->textBrowser->append("2nd Rank = "+QString::number(inter2.get_rank())+"\n");
    QString ker2 = matrix_to_qstring(inter2.get_kernel());
    ui->textBrowser->append("2nd Kernel =");
    ui->textBrowser->append(ker2);

    int main1, main2;
    int gen1, gen2;

    inter1.calcmax();
    inter2.calcmax();

    gen1 = inter1.get_genmax();
    gen2 = inter2.get_genmax();

    if (gen1 == 1)
      main1 = inter1.get_maxa();
    else main1 = inter1.get_maxb();

    if (gen2 == 1)
      main2 = inter2.get_maxa();
    else main2 = inter2.get_maxb();

    if (main1 >= main2) {
      if(gen1 == 1){
        main = 'x'; hidden = 'y';
      }
      else{
        main = 'y'; hidden = 'x';
      }
      f1 = new BivPoly(inter1.generatePoly(gen1));
      f2 = new BivPoly(inter2.generatePoly(gen1));
    }
    else {
      if(gen1 == 1){
        main = 'x'; hidden = 'y';
      }
      else{
        main = 'y'; hidden = 'x';
      }
        f1 = new BivPoly(inter1.generatePoly(gen2));
        f2 = new BivPoly(inter2.generatePoly(gen2));
    }

    //printing equations in the correct txtbox
    ui->equationsTxt1->append(f1->printfunc());
    ui->equationsTxt2->append(f2->printfunc());
}

/*****************************************************************************
******************************************************************************/
//Solving function,solves the given equations system
void Visualization::solve() {

    int B=-1; //If B=-1 we do not solve
    int i;
    int choice;

    double kappa = -1;

    Sylv sylv1(f1,f2);
    B = 7;
    PolyMatrixSy Pmsy(sylv1);
    choice = Pmsy.numofmatr();
    Pmsy.pascal(choice-1);

    // If B==-1 we do not solve the problem
    if( B == -1){
        return;
    }

    kappa = Pmsy.findkappa();
    MatrixXd realxy;
    MatrixXd realmultiys;
    bool flag = true;
    int multiycount = 0;

    int counter;

    if (kappa <= pow(10,B)) {

      ui->textBrowser->append("\nKappa is = "+QString::number(kappa)+" < Bound: "+QString::number(pow(10,B))+" , Standard Eigenproblem\n\n");

      Companion compan(Pmsy);
      realxy = compan.solver();
      counter = realxy.rows();
      if (compan.get_multiy()) {
            realmultiys = compan.get_multiys();
      }
      else {
        flag = false;
        ui->textBrowser->append("No y's with multiplicity > 1 where found\n");
      }
    }
    else {
        ui->textBrowser->append("\nKappa is = "+QString::number(kappa)+" > Bound: "+QString::number(pow(10,B))+" , Generalized Eigenproblem\n\n");

        Generalized generalized(Pmsy);
        realxy = generalized.solver();
        counter = realxy.rows();
        if (generalized.get_multiy()) {
            realmultiys = generalized.get_multiys();
        }
        else {
          flag = false;
          ui->textBrowser->append("No y's with multiplicity > 1 were found\n");
        }
    }

    UniPoly *upol1;
    UniPoly *upol2;
    Companion *comp1;
    Companion *comp2;
    MatrixXd temp1,comptemp1;
    MatrixXd temp2,comptemp2;

    //If y with multiplicity were found then find their x
    if(flag){
        multiycount = realmultiys.rows();
        for (int i = 0; i < multiycount; i++) {
            temp1 = f1->solvey(realmultiys(i));
            temp2 = f2->solvey(realmultiys(i));
            upol1 = new UniPoly(temp1);
            upol2 = new UniPoly(temp2);

            if(f1->get_maxuni() == 0 || f2->get_maxuni() == 0){
                //ui->textBrowser->append("These functions have no roots for y = "+QString::number(realmultiys(i))+"\n");
                continue;
            }
            comp1 = new Companion(upol1,f1->get_maxuni());
            comp2 = new Companion(upol2,f2->get_maxuni());
            comptemp1 = comp1->unisolver();
            comptemp2 = comp2->unisolver();

            bool solcheck1(false);
            bool solcheck2(false);

            for(int k = 0; k < comptemp1.rows(); k++){
              for(int m = 0; m < comptemp2.rows(); m++){
                if( fabs(comptemp1(k) - comptemp2(m)) < pow(10,-5) ){
                   // ui->textBrowser->append("These functions have no roots for y = "+QString::number(realmultiys(i))+"\n");

                    if(main == 'x'){
                        solcheck1 = f1->solveReal(comptemp1(k),realmultiys(i));
                        solcheck2 = f2->solveReal(comptemp1(k),realmultiys(i));
                    }
                    else{
                        solcheck1 = f1->solveReal(realmultiys(i),comptemp1(k));
                        solcheck2 = f2->solveReal(realmultiys(i),comptemp1(k));
                    }

                    if(solcheck1 && solcheck2){
                        //cout << "y = " << realmultiys(i) << " x = " << comptemp1(k) << "    CONFIRMED" << endl;
                        if(abs(comptemp1(k)) > max_x) max_x = abs(comptemp1(k));
                        if(abs(realmultiys(i)) > max_y) max_y = abs(realmultiys(i));
                        ui->textBrowser->append("Root: x = "+QString::number(comptemp1(k))+" y = "+QString::number(realmultiys(i))+"  Confirmed !\n");
                        output << comptemp1(k) << " " << realmultiys(i) << " " << 0 << endl;
                        roots_count++;
                    }
                }
              }
            }

            delete comp1;
            delete comp2;
            delete upol1;
            delete upol2;
      }
    }

    bool solcheck1(false);
    bool solcheck2(false);

    //Solving starting functions with real solutions from eigen
    //ui->textBrowser->append("\nReal roots that have been checked and confirmed are:\n");

    for(i = 0; i < counter; i++){
        //cout << realxy(i,0) << "   " << realxy(i,1) << endl;
        if(main == 'x'){
            solcheck1 = f1->solveReal(realxy(i,1),realxy(i,0));
            solcheck2 = f2->solveReal(realxy(i,1),realxy(i,0));
        }
        else{
            solcheck1 = f1->solveReal(realxy(i,0),realxy(i,1));
            solcheck2 = f2->solveReal(realxy(i,0),realxy(i,1));
        }

        if(solcheck1 && solcheck2){
            //ui->textBrowser->append("y = "+QString::number(realxy(i,0))+"  x = "+QString::number(realxy(i,1))+"\n");
            if(abs(realxy(i,1)) > max_x) max_x = abs(realxy(i,1));
            if(abs(realxy(i,0)) > max_y) max_y = abs(realxy(i,0));
            ui->textBrowser->append("Root: x = "+QString::number(realxy(i,1))+" y = "+QString::number(realxy(i,0))+"  Confirmed !\n");
            output << realxy(i,1) << " " << realxy(i,0) << " " << 0 << endl;
            roots_count++;
        }
    }
    output << endl << endl;

    /**************************** VARIABLE CHANGE ****************************
    *************************************************************************/
    //If variable change option is clicked then follow the algorithm
    //as done in 2nd exercise
    if(varchange){

        ui->textBrowser->append("###### Variable Change Results ######");

        int matrixsize = sylv1.get_total();
        int kappaNew;

        bool newflag = true;
        int newmultiycount = 0;
        int newcounter;
        MatrixXd newrealxy;
        MatrixXd newrealmultiys;

        PolyMatrixSy PmsyNew(choice,matrixsize);
        PmsyNew.pascal(choice-1);

        Pmsy.varchange(PmsyNew);

        // choice = PmsyNew.numofmatr();
        // cout << "New Matrix-Coefficient Md'" << endl;
        // PmsyNew.printindi(choice-1);

        kappaNew = PmsyNew.findkappa();

        // cout << "\nKappa is = " << kappaNew;
        if (kappaNew <= pow(10,B)) {

            ui->textBrowser->append("\nNew Kappa is = "+QString::number(kappaNew)+" < Bound: "+QString::number(pow(10,B))+" , Standard Eigenproblem\n\n");
            // cout << " < Bound: "<<pow(10,B)<< " , Standard Eigenproblem"<<endl<<endl;

            Companion companNew(PmsyNew);
            newrealxy = companNew.solver();
            newcounter = newrealxy.rows();

            if (companNew.get_multiy()) {
                  newrealmultiys = companNew.get_multiys();
                  newmultiycount = newrealmultiys.rows();
            }
            else {
              newflag = false;
              ui->textBrowser->append("No new y's with multiplicity > 1 were found\n");
            }
        }
        else {

            //cout << " > Bound: "<<pow(10,B)<<" , Generalized Eigenproblem"<<endl<<endl;
            ui->textBrowser->append("\nNew Kappa' is = "+QString::number(kappaNew)+" > Bound: "+QString::number(pow(10,B))+" , Generalized Eigenproblem\n\n");
            Generalized generalizedNew(PmsyNew);
            newrealxy = generalizedNew.solver();
            newcounter = newrealxy.rows();

            if (generalizedNew.get_multiy()) {
                newrealmultiys = generalizedNew.get_multiys();
                newmultiycount = newrealmultiys.rows();
            }
            else {
              newflag = false;
              ui->textBrowser->append("No new y's with multiplicity > 1 were found\n");
            }
        }

        if(newflag){
            for (int c = 0; c < newmultiycount; c++){
                //ui->textBrowser->append("New Root: x' = "+QString::number(newrealxy(i,1))+" y' = "+QString::number(newrealxy(i,0))+" Found !\n");
                ui->textBrowser->append("New multiple Root: y' = "+QString::number(newrealmultiys(c))+" Found !\n");
            }
        }

        for(i = 0; i < newcounter; i++){
            ui->textBrowser->append("New Root: x' = "+QString::number(newrealxy(i,1))+" y' = "+QString::number(newrealxy(i,0))+" Found !\n");
            //cout << realxy(i,0) << "   " << realxy(i,1) << endl;
        }
    }
}

/*****************************************************************************
******************************************************************************/
//Printing given points in the correct text widget
void Visualization::printall() {
    if (flagk == 1) {
        ui->listWidget_1->clear();
        ui->listWidget_1->addItems( points1 );
    }
    else
        ui->listWidget_2->clear();
        ui->listWidget_2->addItems( points2 );
}


/******************************************************************************
******                      UI Buttons Functions                         ******
******************************************************************************/
// Button for set & plot 1st function points
void Visualization::on_pushButton_clicked()
{
    QString d = ui->textEdit_1->toPlainText();
    deg1 = d.toInt();
    flagk = 1;
    k1 = (((deg1 + 1) * (deg1 + 2)) / 2) - 1;

    points1.clear();
    widget->show();
    widget->resize(QSize(800,600));

    QString aa = "set tics scale 0.75\nset xtics 1\nset ytics 1\nset xrange [";
    QString bb = "]\nset yrange [";
    QString cc = "]\nset xlabel 'x'\nset ylabel 'y'\nset zeroaxis\nplot \"<echo \" notitle\n";
    QString final = aa+xrange+bb+yrange+cc;

    *instance << final;
}

/*****************************************************************************
******************************************************************************/
// Button for set & plot 2nd function points
void Visualization::on_pushButton_3_clicked()
{
    QString d = ui->textEdit_2->toPlainText();
    deg2 = d.toInt();
    flagk = 2;
    k2 = (((deg2 + 1) * (deg2 + 2)) / 2) - 1;

    points2.clear();
    widget->show();
    widget->resize(QSize(800,600));

    QString aa = "set tics scale 0.75\nset xtics 1\nset ytics 1\nset xrange [";
    QString bb = "]\nset yrange [";
    QString cc = "]\nset xlabel 'x'\nset ylabel 'y'\nset zeroaxis\nplot \"<echo \" notitle\n";
    QString final = aa+xrange+bb+yrange+cc;

    *instance << final;
}

/*****************************************************************************
******************************************************************************/
//Generate equations button from Interpolation method tab
void Visualization::on_pushButton_6_clicked()
{
    generate_interpolation();
}

/*****************************************************************************
******************************************************************************/
//Solve button which calls the solving function
void Visualization::on_pushButton_7_clicked()
{
    solve();
}

/*****************************************************************************
******************************************************************************/
//Plot button which shows the final plot
void Visualization::on_pushButton_8_clicked()
{
    QString gnuf1 = f1->printfunc().replace("^", "**");
    gnuf1.replace(" ", "");
    QString gnuf2 = f2->printfunc().replace("^", "**");
    gnuf2.replace(" ", "");

    for (int i = 0; i < pts1.rows(); i++) {
        output << pts1(i,0) << " " << pts1(i,1) << " " << 0 << endl;
    }
    output << endl << endl;
    for (int i = 0; i < pts2.rows(); i++) {
        output << pts2(i,0) << " " << pts2(i,1) << " " << 0 << endl;
    }

    if(roots_count != 0){
        max_x = ceil(max_x) + 5;
        max_y = ceil(max_y) + 5;
    }

    sol_xrange = ("-"+QString::number(max_x)+":"+QString::number(max_x));
    sol_yrange = ("-"+QString::number(max_y)+":"+QString::number(max_y));

    widget->show();
    widget->resize(QSize(800,600));
    QString pp;
    QString aa = "set yrange [";
    QString nn = "]\nset xrange [";
    QString mm = "]\nset isosamples 500,500\nf1(x,y)= ";
    QString s3 ="\nf2(x,y)= ";

    if(roots_count == 0){
        pp = "\nset contour\nset cntrparam levels discrete 0\nset view 0,0\nunset ztics\nset surface\nset style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 0.75"
             "\nset style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 0.75\nsplot 'roots' index 0 with points ls 1 title '1st eq. points', \''index 1 with points ls 2 title '2nd eq. points' nocontour, f1(x,y) with lines nosurface ,f2(x,y) with lines nosurface  \n";
    }
    else{
        pp = "\nset contour\nset cntrparam levels discrete 0\nset view 0,0\nunset ztics\nset surface\nset style line 1 lc rgb '#ff9500' lt 1 lw 2 pt 1 ps 1.5\nset style line 2 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 0.75"
             "\nset style line 3 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 0.75\nsplot 'roots' index 0 with points ls 1 title 'Common Roots', \''index 1 with points ls 2 title '1st eq. points', \''index 2 with points ls 3 title '2nd eq. points' nocontour, f1(x,y) with lines nosurface ,f2(x,y) with lines nosurface  \n";
    }

    QString cc = aa+sol_yrange+nn+sol_xrange+mm+gnuf1+s3+gnuf2+pp;
    *instance << cc;
}

/*****************************************************************************
******************************************************************************/
//Range input set button for points plotting
void Visualization::on_pushButton_9_clicked()
{
    xrange = ui->textEdit_3->toPlainText();
    max_x = xrange.toInt();
    xrange.prepend("-"+xrange+":");

    yrange = ui->textEdit_4->toPlainText();
    max_y = yrange.toInt();
    yrange.prepend("-"+yrange+":");
}

/*****************************************************************************
******************************************************************************/
//File open(browse) button from Read File Tab
void Visualization::on_pushButton_2_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QString(),
                tr("All Files (*.*)"));

    ui->textBrowser_2->setPlainText(fileName);
    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly)) {
            QMessageBox::critical(this, tr("Error"), tr("Could not open file"));
            return;
        }
        QTextStream in(&file);
        QString str1 = in.readLine();
        QString str2 = in.readLine();
        func1 = str1.toUtf8().constData();
        func2 = str2.toUtf8().constData();
        ui->equationsTxt1->setText(str1);
        ui->equationsTxt2->setText(str2);
        file_read();
        file.close();
    }
}

/*****************************************************************************
******************************************************************************/
//Checkbox for variable change option
void Visualization::on_checkBox_toggled(bool checked)
{
    if( checked ){
        varchange = true;
        //ui->textBrowser->append(" JUST CHECKED ");
    }
    else{
        varchange = false;
        //ui->textBrowser->append(" JUST UNCHECKED ");
    }
}

/*****************************************************************************
******************************************************************************/
// Deg1 from Read File Tab
void Visualization::on_pushButton_4_clicked()
{
    QString d = ui->textEdit_7->toPlainText();
    deg1 = d.toInt();
}

/*****************************************************************************
******************************************************************************/
// Deg2 from Read File Tab
void Visualization::on_pushButton_5_clicked()
{
    QString d = ui->textEdit_6->toPlainText();
    deg2 = d.toInt();
}


/*****************************************************************************
 *******                    Help Functions                            ********
******************************************************************************/

//Converts the given QStringList to matrix
MatrixXd qlist_to_matrix(QStringList lst) {
    QString strng = lst.join(",");
    QStringList strlst = strng.split(",");
    MatrixXd temp(lst.size(),2);
    int j = 0;
    for (int i = 0; i < strlst.size(); i += 2) {
        temp(j,0) = strlst.at(i).toDouble();
        temp(j,1) = strlst.at(i+1).toDouble();
        j++;
    }
    return temp;
}

/*****************************************************************************
******************************************************************************/
//Converts the given matrix to QString
QString matrix_to_qstring(MatrixXd mtrx) {
    QString str;
    for (int i = 0; i < mtrx.rows(); i++) {
        for (int j = 0; j < mtrx.cols(); j++) {
            str.append(QString::number(mtrx(i,j)));
            str.append("  ");
        }
        str.append("\n");
    }
    return str;
}
