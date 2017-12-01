#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#ifdef HAVE_QT5
#include <QtWidgets/QMainWindow>
#else
#include <QtGui/QMainWindow>
#endif

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
	void inputFileSelecter();
    void siteFileSelecter();
    void outputFileSelecter();
    void populateTable();
    void processLine(QString line);
//    void join_checkbox();
//    void savepaths_checkbox();
    void on_radioButton_len_limit_toggled(bool checked);
    void on_radioButton_len_diam_limit_toggled(bool checked);
    void on_checkBoxImage_toggled(bool checked);
    void batch_analyser();
    int blocker(QString numstr);
    int distancer(QString numstr);
    int tiffer(QString numstr);
    int reader(QString ca_outputfile);

private:
    Ui::MainWindow *ui;

public:
	QString inputFileName;
    QString outputBaseFileName;
    QString siteFileName;
    QString subregionFileName;
    QString tempfileName;
    QStringList headers;
    QStringList variableNames;
    QStringList calcNames;
    QStringList values;
    QString x0_str, y0_str, z0_str, R_str;
};

#endif // MAINWINDOW_H
