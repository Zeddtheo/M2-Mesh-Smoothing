#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>

//-------------------------------TP-Lissage Laplacien de Maillages------------------------------------
// Calcul du cotangent d'un angle
float MainWindow::cot(float angle){
    return 1.0f/(std::tan(angle));
}

// Calcul de l'aire d'un triangle à partir de ses sommets
double MainWindow::calc_area(MyMesh::Point p[]){
    MyMesh::Point v1 = p[0];
    MyMesh::Point v2 = p[1];
    MyMesh::Point v3 = p[2];

    MyMesh::Point u = v2 - v1;
    MyMesh::Point v = v3 - v1;

    double ux = u[0], uy = u[1], uz = u[2];
    double vx = v[0], vy = v[1], vz = v[2];

    double i = (uy*vz - vy*uz);
    double j = (ux*vz - vx*uz);
    double k = (ux*vy - vx*uy);

    double aire = i - j + k;
    double area = abs(aire)/2.0;
    return area;
}

// Compte le nombre d'arêtes d'un sommet dans un maillage
int MainWindow::count_vertex_edges(MyMesh* _mesh, MyMesh::VertexHandle vh) {
    int count = 0;

    for (MyMesh::VertexEdgeIter ve_it = _mesh->ve_iter(vh); ve_it.is_valid(); ++ve_it) {
        ++count;
    }
    return count;
}

// Calcule la somme des aires des faces voisines d'un sommet
double MainWindow::calc_neighbour_area(MyMesh* _mesh, VertexHandle v){
    double somme = 0;
    double aire = 0;
    MyMesh::Point points[3];
    int cpt= 0;
    for(MyMesh::VertexFaceIter f = _mesh->vf_iter(v); f.is_valid(); f++){
        cpt = 0;
        for(MyMesh::FaceVertexIter fv_it = _mesh->fv_iter(*f); fv_it.is_valid(); fv_it++){
            points[cpt++] = _mesh->point(*fv_it);
        }
        aire = calc_area(points);
        somme += aire;
    }

    return somme/3;
}

// Calcul des coefficients de cotangente pour un sommet
Vec3f MainWindow::calc_cot_coef(MyMesh *_mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    Vec3f f_v = Vec3f(0.0f, 0.0f, 0.0f);
    Vec3f v_pos = _mesh->point(vh);

    for (MyMesh::VertexOHalfedgeIter voh_it = _mesh->voh_iter(vh); voh_it.is_valid(); ++voh_it){
        HalfedgeHandle v_vi = *voh_it;
        MyMesh::VertexHandle vi = _mesh->to_vertex_handle(v_vi);
        MyMesh::Point vi_pos = _mesh->point(vi);

        MyMesh::HalfedgeHandle vi_b = _mesh->next_halfedge_handle(v_vi);
        MyMesh::VertexHandle b = _mesh->to_vertex_handle(vi_b);

        MyMesh::HalfedgeHandle vi_v = _mesh->opposite_halfedge_handle(v_vi);
        MyMesh::HalfedgeHandle v_a = _mesh->next_halfedge_handle(vi_v);
        MyMesh::VertexHandle a = _mesh->to_vertex_handle(v_a);

        MyMesh::Point v3 = _mesh->point(vh);
        MyMesh::Point vi3 = _mesh->point(vi);
        MyMesh::Point ai3 = _mesh->point(a);
        MyMesh::Point bi3 = _mesh->point(b);
        // On utilise les halfEdge pour trouver les 4 vecteurs qui nous interessent
        //pour pouvoir calculer les angles alpha et beta
        Vec3f a_v = v3 - ai3;
        Vec3f a_vi = vi3 - ai3;
        Vec3f b_v = v3 - bi3;
        Vec3f b_vi = vi3 - bi3;

        float angle_a = acos(OpenMesh::dot(a_v.normalized(), a_vi.normalized()));
        float angle_b = acos(OpenMesh::dot(b_v.normalized(), b_vi.normalized()));

        float cot_angle_a = cot(angle_a);
        float cot_angle_b = cot(angle_b);
        //std::cout<<"somme = "<<cot_angle_a + cot_angle_b<<std::endl;
        f_v += (cot_angle_a + cot_angle_b) * (vi_pos - v_pos);

    }
    //std::cout<<"f_v cot: "<<f_v<<std::endl;
    return f_v;
}

// Calcul des coefficients uniformes pour un sommet
Vec3f MainWindow::calc_coef_uniforme(MyMesh *_mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    Vec3f f_v(0.0f, 0.0f, 0.0f);
    Vec3f v_pos = _mesh->point(vh);

    for (MyMesh::VertexOHalfedgeIter voh_it = _mesh->voh_iter(vh); voh_it.is_valid(); ++voh_it){
        HalfedgeHandle v_vi = *voh_it;
        MyMesh::VertexHandle vi_vh = _mesh->to_vertex_handle(v_vi);
        Vec3f vi_pos = _mesh->point(vi_vh);
        f_v += vi_pos - v_pos;
        //std::cout<<"uniforme: "<<(vi_pos - v_pos)*100<<std::endl;
    }
    //std::cout<<"f_v uniforme: "<<f_v<<std::endl;
    return f_v;
}

// Approximation cotangent pour le lissage de Laplace
std::vector<Vec3f> MainWindow::approximation_cotgentielle(MyMesh* _mesh) {
    std::vector<Vec3f> operator_laplace(_mesh->n_vertices());

    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!= mesh.vertices_end(); ++v_it){
        //double aire = calc_area_barycentric(_mesh, v_it->idx());
        MyMesh::VertexHandle vh = *v_it;
        double aire = calc_neighbour_area(_mesh, vh);
        //std::cout<<"aire: "<<aire<<std::endl;
        //float area = 1/(2.0*calc_area_barycentric(_mesh, v_it->idx()));

        MyMesh::Point weight = calc_cot_coef(_mesh, v_it->idx());
        //std::cout<<"weight: "<<weight<<std::endl;
        MyMesh::Point result = 1/(2.0 * aire) * weight;

        operator_laplace.at(v_it->idx()) = result;
        operator_laplace.at(v_it->idx()).normalize();

        //std::cout<<operator_laplace.at(v_it->idx())<<std::endl;
    }
    return operator_laplace;
}

// Approximation uniforme pour le lissage de Laplace
std::vector<Vec3f> MainWindow::approximation_uniforme(MyMesh* _mesh) {
    std::vector<Vec3f> operator_laplace(_mesh->n_vertices());

    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!= mesh.vertices_end(); ++v_it){
        MyMesh::VertexHandle vh = *v_it;
        Vec3f weight = calc_coef_uniforme(_mesh, v_it->idx());
        Vec3f result = (1.0f / count_vertex_edges(_mesh, vh)) * weight;
        operator_laplace.at(v_it->idx()) = result;
    }

    return operator_laplace;
}

// Calcul du poids dans la matrice pour les sommets voisins
double MainWindow::calc_w(MyMesh *_mesh, VertexHandle v, VertexHandle vi){
    HalfedgeHandle half_edge;
    for (MyMesh::VertexOHalfedgeIter he = _mesh->voh_iter(v); he.is_valid(); he++) {
        if (_mesh->to_vertex_handle(*he) == vi) {
            half_edge = *he;
            break;
        }
    }

    MyMesh::Point p_v = _mesh->point(v);
    MyMesh::Point p_vi = _mesh->point(vi);

    HalfedgeHandle next = _mesh->next_halfedge_handle(half_edge);
    VertexHandle vh_alpha = _mesh->to_vertex_handle(next);
    MyMesh::Point p_alpha = _mesh->point(vh_alpha);
    Vec3f a_v = p_v - p_alpha;
    Vec3f a_vi = p_vi - p_alpha;
    float angle_alpha = acos(OpenMesh::dot(a_v.normalized(), a_vi.normalized()));

    half_edge = _mesh->opposite_halfedge_handle(half_edge);
    next = _mesh->next_halfedge_handle(half_edge);
    VertexHandle vh_beta = _mesh->to_vertex_handle(next);
    MyMesh::Point p_beta = _mesh->point(vh_beta);
    Vec3f b_v = p_v - p_beta;
    Vec3f b_vi = p_vi - p_beta;
    float angle_beta = acos(OpenMesh::dot(b_v.normalized(), b_vi.normalized()));

    return cot(angle_alpha) + cot(angle_beta);
}

//Matrix M
Eigen::SparseMatrix<double> MainWindow::matrix_D(MyMesh* _mesh) {
    int n = _mesh->n_vertices();
    Eigen::SparseMatrix<double> d(n, n);

    for (MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); ++v) {
        VertexHandle vh = *v;
        double aire = calc_neighbour_area(_mesh, vh);
        d.insert(v->idx(), v->idx()) = 1.0 / (2.0 * aire);
    }

    return d;
}

// Matrix D
Eigen::SparseMatrix<double> MainWindow::matrix_M(MyMesh* _mesh) {
    int nb_sommets = _mesh->n_vertices();
    Eigen::SparseMatrix<double> m(nb_sommets, nb_sommets);

    for (MyMesh::VertexIter v_i = _mesh->vertices_begin(); v_i != _mesh->vertices_end(); ++v_i) {
        VertexHandle v = *v_i;
        double weight = 0.0;

        for (MyMesh::VertexVertexIter vv_it = _mesh->vv_iter(v); vv_it.is_valid(); ++vv_it) {
            VertexHandle vi = *vv_it;

            if (v != vi) {
                double w = calc_w(_mesh, v, vi);
                m.insert(v.idx(), vi.idx()) = w;
                weight += w;
            }
        }

        m.insert(v.idx(), v.idx()) = -weight;
    }
    return m;
}

// Matrix L
Eigen::SparseMatrix<double> MainWindow::matrix_L(MyMesh* _mesh) {

    Eigen::SparseMatrix<double> m = matrix_M(_mesh);
    Eigen::SparseMatrix<double> d = matrix_D(_mesh);

    Eigen::SparseMatrix<double> l_matrix = d * m;

    return l_matrix;
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

void MainWindow::on_pushButton_lissage_clicked() {
    float h = 0.1f;
    float lambda = 0.03f;
    std::vector<Vec3f> mesh_lissage = approximation_cotgentielle(&mesh);

    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        mesh.point(mesh.vertex_handle(v_it->idx())) += h * lambda * mesh_lissage[v_it->idx()];
    }
    displayMesh(&mesh);
}

void MainWindow::on_pushButton_lissage_uniforme_clicked() {
    float h_uniforme = 0.5f;
    float l_uniforme = 0.5f;
    std::vector<Vec3f> mesh_lissage = approximation_uniforme(&mesh);


    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        mesh.point(mesh.vertex_handle(v_it->idx())) += h_uniforme * l_uniforme * mesh_lissage[v_it->idx()];
    }
    displayMesh(&mesh);
}

void MainWindow::on_pushButton_matrix_clicked() {
    Eigen::SparseMatrix<double> lb_matrix = matrix_L(&mesh);

    std::cout << "Laplace-Beltrami Matrix:" << std::endl;
    for (int i = 0; i < lb_matrix.rows(); ++i) {
        for (int j = 0; j < lb_matrix.cols(); ++j) {
            if (lb_matrix.coeff(i, j) != 0) {
                std::cout << "lb_matrix(" << i << ", " << j << ") = " << lb_matrix.coeff(i, j) << std::endl;
            }
        }
    }
}

//-------------------------------Fin de TP------------------------------------
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh *_mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh *_mesh, DisplayMode mode)
{
    GLuint *triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat *triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat *triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if (mode == DisplayMode::TemperatureMap)
    {
        QVector<float> values;
        for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
            values.append(fabs(_mesh->data(*curVert).value));
        std::sort(values.begin(), values.end());

        float range = values.at(values.size() * 0.8);

        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt != fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if (_mesh->data(*fvIt).value > 0)
            {
                triCols[3 * i + 0] = 255;
                triCols[3 * i + 1] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 2] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            else
            {
                triCols[3 * i + 2] = 255;
                triCols[3 * i + 1] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 0] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            if (_mesh->data(*fvIt).value > 0)
            {
                triCols[3 * i + 0] = 255;
                triCols[3 * i + 1] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 2] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            else
            {
                triCols[3 * i + 2] = 255;
                triCols[3 * i + 1] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 0] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            if (_mesh->data(*fvIt).value > 0)
            {
                triCols[3 * i + 0] = 255;
                triCols[3 * i + 1] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 2] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            else
            {
                triCols[3 * i + 2] = 255;
                triCols[3 * i + 1] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 0] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if (mode == DisplayMode::Normal)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt != fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3 * i + 0] = _mesh->color(*fIt)[0];
            triCols[3 * i + 1] = _mesh->color(*fIt)[1];
            triCols[3 * i + 2] = _mesh->color(*fIt)[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            triCols[3 * i + 0] = _mesh->color(*fIt)[0];
            triCols[3 * i + 1] = _mesh->color(*fIt)[1];
            triCols[3 * i + 2] = _mesh->color(*fIt)[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            triCols[3 * i + 0] = _mesh->color(*fIt)[0];
            triCols[3 * i + 1] = _mesh->color(*fIt)[1];
            triCols[3 * i + 2] = _mesh->color(*fIt)[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if (mode == DisplayMode::ColorShading)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt != fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3 * i + 0] = _mesh->data(*fvIt).faceShadingColor[0];
            triCols[3 * i + 1] = _mesh->data(*fvIt).faceShadingColor[1];
            triCols[3 * i + 2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            triCols[3 * i + 0] = _mesh->data(*fvIt).faceShadingColor[0];
            triCols[3 * i + 1] = _mesh->data(*fvIt).faceShadingColor[1];
            triCols[3 * i + 2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            triCols[3 * i + 0] = _mesh->data(*fvIt).faceShadingColor[0];
            triCols[3 * i + 1] = _mesh->data(*fvIt).faceShadingColor[1];
            triCols[3 * i + 2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint *linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat *linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat *linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int>> edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if (t > 0)
        {
            if (!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int>> it(edgesIDbyThickness);
    QList<QPair<float, int>> edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for (int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3 * i + 0] = _mesh->point(vh1)[0];
            linesVerts[3 * i + 1] = _mesh->point(vh1)[1];
            linesVerts[3 * i + 2] = _mesh->point(vh1)[2];
            linesCols[3 * i + 0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3 * i + 1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3 * i + 2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3 * i + 0] = _mesh->point(vh2)[0];
            linesVerts[3 * i + 1] = _mesh->point(vh2)[1];
            linesVerts[3 * i + 2] = _mesh->point(vh2)[2];
            linesCols[3 * i + 0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3 * i + 1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3 * i + 2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint *pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat *pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat *pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int>> vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if (t > 0)
        {
            if (!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int>> vitt(vertsIDbyThickness);
    QList<QPair<float, int>> vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for (int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3 * i + 0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3 * i + 1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3 * i + 2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3 * i + 0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3 * i + 1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3 * i + 2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


