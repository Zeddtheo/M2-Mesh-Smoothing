/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 6.5.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <meshviewerwidget.h>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QWidget *widget_2;
    QGridLayout *gridLayout;
    QVBoxLayout *verticalLayout;
    QPushButton *pushButton_chargement;
    QPushButton *pushButton_lissage;
    QPushButton *pushButton_lissage_uniforme;
    QPushButton *pushButton_matrix;
    QPushButton *pushButton_H;
    QPushButton *pushButton_K;
    QPushButton *pushButton_angleArea;
    MeshViewerWidget *displayWidget;
    QMenuBar *menuBar;
    QMenu *menuEditer;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName("MainWindow");
        MainWindow->resize(632, 408);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName("centralWidget");
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName("horizontalLayout");
        widget_2 = new QWidget(centralWidget);
        widget_2->setObjectName("widget_2");
        QSizePolicy sizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(150);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(widget_2->sizePolicy().hasHeightForWidth());
        widget_2->setSizePolicy(sizePolicy);
        widget_2->setMinimumSize(QSize(150, 0));
        gridLayout = new QGridLayout(widget_2);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName("gridLayout");
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName("verticalLayout");
        pushButton_chargement = new QPushButton(widget_2);
        pushButton_chargement->setObjectName("pushButton_chargement");
        pushButton_chargement->setMinimumSize(QSize(200, 0));

        verticalLayout->addWidget(pushButton_chargement);

        pushButton_lissage = new QPushButton(widget_2);
        pushButton_lissage->setObjectName("pushButton_lissage");

        verticalLayout->addWidget(pushButton_lissage);

        pushButton_lissage_uniforme = new QPushButton(widget_2);
        pushButton_lissage_uniforme->setObjectName("pushButton_lissage_uniforme");

        verticalLayout->addWidget(pushButton_lissage_uniforme);

        pushButton_matrix = new QPushButton(widget_2);
        pushButton_matrix->setObjectName("pushButton_matrix");

        verticalLayout->addWidget(pushButton_matrix);

        pushButton_H = new QPushButton(widget_2);
        pushButton_H->setObjectName("pushButton_H");

        verticalLayout->addWidget(pushButton_H);

        pushButton_K = new QPushButton(widget_2);
        pushButton_K->setObjectName("pushButton_K");

        verticalLayout->addWidget(pushButton_K);

        pushButton_angleArea = new QPushButton(widget_2);
        pushButton_angleArea->setObjectName("pushButton_angleArea");

        verticalLayout->addWidget(pushButton_angleArea);


        gridLayout->addLayout(verticalLayout, 0, 0, 1, 1);


        horizontalLayout->addWidget(widget_2);

        displayWidget = new MeshViewerWidget(centralWidget);
        displayWidget->setObjectName("displayWidget");
        displayWidget->setAutoFillBackground(false);

        horizontalLayout->addWidget(displayWidget);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName("menuBar");
        menuBar->setGeometry(QRect(0, 0, 632, 22));
        menuEditer = new QMenu(menuBar);
        menuEditer->setObjectName("menuEditer");
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName("mainToolBar");
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName("statusBar");
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuEditer->menuAction());

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        pushButton_chargement->setText(QCoreApplication::translate("MainWindow", "Charger OBJ", nullptr));
        pushButton_lissage->setText(QCoreApplication::translate("MainWindow", "Lissage", nullptr));
        pushButton_lissage_uniforme->setText(QCoreApplication::translate("MainWindow", "Lissage uniforme", nullptr));
        pushButton_matrix->setText(QCoreApplication::translate("MainWindow", "ShowMatrix", nullptr));
        pushButton_H->setText(QCoreApplication::translate("MainWindow", "Courbure Moyenne", nullptr));
        pushButton_K->setText(QCoreApplication::translate("MainWindow", "Courbure Gaussienne", nullptr));
        pushButton_angleArea->setText(QCoreApplication::translate("MainWindow", "Test angles/aires", nullptr));
        menuEditer->setTitle(QCoreApplication::translate("MainWindow", "Editer", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
