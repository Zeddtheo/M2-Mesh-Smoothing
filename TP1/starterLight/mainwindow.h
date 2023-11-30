#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness; float value; Color faceShadingColor;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

enum DisplayMode {Normal, TemperatureMap, ColorShading};

class MainWindow : public QMainWindow
{
    Q_OBJECT
   // ui->setupUi(this);
//    connect(ui->sliderH, &QSlider::valueChanged, this, &MyWindow::updateH);
//    connect(ui->sliderLambda, &QSlider::valueChanged, this, &MyWindow::updateLambda);
    //connect(ui->doubleSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::onDoubleSpinBoxValueChanged);

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    //Fonctions lissage
    float cot(float angle);
    int count_vertex_edges(MyMesh* _mesh, MyMesh::VertexHandle vhandle);

    double calc_neighbour_area(MyMesh* _mesh, VertexHandle v);
    double calc_area(MyMesh::Point p[]);

    Vec3f calc_cot_coef(MyMesh *_mesh, int vertexID);
    Vec3f calc_coef_uniforme(MyMesh *_mesh, int vertexID);

    std::vector<Vec3f> approximation_cotgentielle(MyMesh* _mesh);
    std::vector<Vec3f> approximation_uniforme(MyMesh* _mesh);

    //Fonctions matrix calcul
    double calc_w(MyMesh *_mesh, VertexHandle v, VertexHandle vi);
    Eigen::SparseMatrix<double> matrix_D(MyMesh* _mesh);
    Eigen::SparseMatrix<double> matrix_M(MyMesh* _mesh);
    Eigen::SparseMatrix<double> matrix_L(MyMesh* _mesh);

    //Fonctions mesh affichage
    void displayMesh(MyMesh *_mesh, DisplayMode mode = DisplayMode::Normal);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:
    void on_pushButton_chargement_clicked();
    void on_pushButton_lissage_clicked();
    void on_pushButton_lissage_uniforme_clicked();
    void on_pushButton_matrix_clicked();

private:
    bool modevoisinage;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
