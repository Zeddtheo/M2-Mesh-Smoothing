#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#define PI 3.14159265


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
//    ui->setupUi(this);
//    connect(ui->sliderH, &QSlider::valueChanged, this, &MyWindow::updateH);
//    connect(ui->sliderLambda, &QSlider::valueChanged, this, &MyWindow::updateLambda);

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();


    // les fonctions à compléter
    float faceArea(MyMesh* _mesh, int faceID);
    float angleFF(MyMesh *_mesh, int faceID0, int faceID1, int vertID0, int vertID1);
    float angleEE(MyMesh* _mesh, int vertexID, int faceID);
    void H_Curv(MyMesh* _mesh);
    void K_Curv(MyMesh* _mesh);
    //TP M2
    float cotan(float angle);
    float calc_area_barycentric(MyMesh* _mesh, int vertexID);
    int count_vertex_edges(MyMesh* _mesh, MyMesh::VertexHandle vhandle);

    Vec3f calc_cot_coef(MyMesh *_mesh, int vertexID);
    Vec3f calc_coef_uniforme(MyMesh *_mesh, int vertexID);

    std::vector<std::vector<float>> LB_Matrix(MyMesh* _mesh);
    std::vector<Vec3f> approximation_cotangentielle(MyMesh* _mesh);
    std::vector<Vec3f> approximation_cotangentielle_uniforme(MyMesh* _mesh);

    bool isNeighbor(int vertexA, int vertexB);

    float calc_sum_angleEE(MyMesh* _mesh, int vertexID);
    float calc_length_edge(MyMesh *_mesh, int vertexID);
    float calc_sum_angleFF(MyMesh *_mesh, int vertexID);

    void displayMesh(MyMesh *_mesh, DisplayMode mode = DisplayMode::Normal);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_angleArea_clicked();
    void on_pushButton_H_clicked();
    void on_pushButton_K_clicked();
    void on_pushButton_lissage_clicked();
    void on_pushButton_lissage_uniforme_clicked();
    void on_pushButton_matrix_clicked();

private:

    bool modevoisinage;
    float h = 0.01f;
    float lambda = 0.02f;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
