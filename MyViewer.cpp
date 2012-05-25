#include <algorithm>
#include <iostream>
#include <vector>

#include <QKeyEvent>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

#include "MyViewer.h"

MyViewer::MyViewer(QWidget *parent) :
  QGLViewer(parent),
  mean_min(0.0), mean_max(0.0), cutoff_ratio(0.05),
  show_mean(false), show_solid(true), show_wireframe(false)
{
  setSelectRegionWidth(5);
  setSelectRegionHeight(5);
}

MyViewer::~MyViewer()
{
}

void MyViewer::updateMeanMinMax()
{
  size_t n = mesh.n_vertices();
  if(n == 0)
    return;

  std::vector<double> mean;
  mean.reserve(n);
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i)  
    mean.push_back(mesh.data(i).mean);

  std::sort(mean.begin(), mean.end());
  size_t k = (double)n * cutoff_ratio;
  mean_min = std::min(mean[k-1], 0.0);
  mean_max = std::max(mean[n-k], 0.0);
}

void MyViewer::updateMeanCurvature(bool update_min_max)
{
  for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i)
    mesh.data(i).area = -1;

  // Compute triangle strip areas
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    mesh.data(i).area = 0;
    mesh.data(i).mean = 0;
    for(MyMesh::ConstVertexFaceIter j(mesh, i); (bool)j; ++j) {
      if(mesh.data(j).area == -1) {
        MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle((OpenMesh::FaceHandle const &)j);
        MyMesh::HalfedgeHandle h2 = mesh.next_halfedge_handle(h1);
        mesh.data(j).area = (halfedgeVector(h1) % halfedgeVector(h2)).norm() / 2.0;
      }
      mesh.data(i).area += mesh.data(j).area;
    }
  }

  // Compute mean values using normal difference angles
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    for(MyMesh::ConstVertexEdgeIter j(mesh, i); (bool)j; ++j) {
      double angle;
      MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle(j, 0);
      MyMesh::HalfedgeHandle h2 = mesh.halfedge_handle(j, 1);
      Vector v = halfedgeVector(h1);
      if(mesh.is_boundary(h1) || mesh.is_boundary(h2))
        angle = 0.0;
      else {
        Vector n1 = mesh.normal(mesh.face_handle(h1));
        Vector n2 = mesh.normal(mesh.face_handle(h2));
        angle = acos(std::min(std::max(n1 | n2, -1.0f), 1.0f));
        angle *= ((n1 % n2) | v) >= 0.0 ? 1.0 : -1.0;
      }
      mesh.data(i).mean += angle * v.norm();
    }
    mesh.data(i).mean *= 3.0 / 4.0 / mesh.data(i).area;
  }

  if(update_min_max)
    updateMeanMinMax();
}

void MyViewer::meanMapColor(double d, double *color) const
{
  if(d <= mean_min) {
    color[0] = 0.0;
    color[1] = 0.0;
    color[2] = 1.0;
  } else if(d >= mean_max) {
    color[0] = 1.0;
    color[1] = 0.0;
    color[2] = 0.0;
  } else if(d < 0) {
    double alpha = d / mean_min;
    color[0] = 0.0;
    color[1] = 1.0 - alpha;
    color[2] = alpha;
  } else {
    double alpha = d / mean_max;
    color[0] = alpha;
    color[1] = 1.0 - alpha;
    color[2] = 0;
  }
}

void MyViewer::fairMesh()
{
  emit startComputation(tr("Fairing mesh..."));
  OpenMesh::Smoother::JacobiLaplaceSmootherT<MyMesh> smoother(mesh);
  smoother.initialize(OpenMesh::Smoother::SmootherT<MyMesh>::Normal, // or: Tangential_and_Normal
                      OpenMesh::Smoother::SmootherT<MyMesh>::C1);
  for(size_t i = 1; i <= 10; ++i) {
    smoother.smooth(10);
    emit midComputation(i * 10);
  }
  mesh.update_face_normals();
  updateMeanCurvature(false);
  emit endComputation();
}

bool MyViewer::openMesh(std::string const &filename)
{
  if(!OpenMesh::IO::read_mesh(mesh, filename) || mesh.n_vertices() == 0)
    return false;
  mesh.request_face_normals();
  mesh.update_face_normals();

  updateMeanCurvature();

  // Set camera on the model
  MyMesh::Point box_min, box_max;
  box_min = box_max = mesh.point(mesh.vertices_begin());
  for(MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    box_min.minimize(mesh.point(i));
    box_max.maximize(mesh.point(i));
  }
  camera()->setSceneBoundingBox(qglviewer::Vec(box_min[0], box_min[1], box_min[2]),
                                qglviewer::Vec(box_max[0], box_max[1], box_max[2]));
  camera()->showEntireScene();

  setSelectedName(-1);
  updateGL();
  return true;
}

void MyViewer::init()
{
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
}

void MyViewer::draw()
{
  if(!show_solid && show_wireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1, 1);

  std::vector<double> color(3, 1.0);
  if(show_solid || show_wireframe)
    for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i) {
      glBegin(GL_POLYGON);
      glNormal3fv(mesh.normal(i).data());
      for(MyMesh::ConstFaceVertexIter j(mesh, i); (bool)j; ++j) {
        if(show_mean)
          meanMapColor(mesh.data(j).mean, &color[0]);
        glColor3dv(&color[0]);
        glVertex3fv(mesh.point(j).data());
      }
      glEnd();
    }

  if(show_solid && show_wireframe) {
    glPolygonMode(GL_FRONT, GL_LINE);
    glColor3d(0.0, 0.0, 0.0);
    glDisable(GL_LIGHTING);
    for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i) {
      glBegin(GL_POLYGON);
      for(MyMesh::ConstFaceVertexIter j(mesh, i); (bool)j; ++j)
        glVertex3fv(mesh.point(j).data());
      glEnd();
    }
    glEnable(GL_LIGHTING);
  }

  if(show_wireframe && selectedName() != -1) {
    glDisable(GL_LIGHTING);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(10.0);
    glColor3d(0.4, 0.4, 1.0);
    glBegin(GL_POINTS);
    glVertex3fv(mesh.point(selected).data());
    glEnd();
    glEnable(GL_LIGHTING);
  }
}

void MyViewer::drawWithNames()
{
  if(!show_wireframe)
    return;

  int j = 0;
  for(MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    glPushName(j++);
    glRasterPos3fv(mesh.point(i).data());
    glPopName();
  }
}

void MyViewer::postSelection(const QPoint &)
{
  int sel = selectedName();
  if(sel == -1)
    return;

  MyMesh::ConstVertexIter i = mesh.vertices_begin();
  for(int j = 0; j != sel; ++i, ++j);
  selected = i;
}

void MyViewer::keyPressEvent(QKeyEvent *e)
{
  if(e->modifiers() == Qt::NoModifier)
    switch(e->key()) {
    case Qt::Key_M:
      show_mean = !show_mean;
      updateGL();
      break;
    case Qt::Key_S:
      show_solid = !show_solid;
      updateGL();
      break;
    case Qt::Key_W:
      show_wireframe = !show_wireframe;
      updateGL();
      break;
    case Qt::Key_F:
      fairMesh();
      updateGL();
      break;
    default:
      QGLViewer::keyPressEvent(e);
    }
  else
    QGLViewer::keyPressEvent(e);
}

void MyViewer::mouseMoveEvent(QMouseEvent *e)
{
  if(selectedName() != -1 && e->modifiers() & Qt::ShiftModifier && e->buttons() & Qt::LeftButton) {
    bool found;
    qglviewer::Vec p = camera()->pointUnderPixel(e->pos(), found);
    if(found) {
      mesh.set_point(selected, MyMesh::Point(p[0], p[1], p[2]));
      updateGL();
    }
  } else
    QGLViewer::mouseMoveEvent(e);
}

QString MyViewer::helpString() const
{
  QString text("<h2>Sample Framework</h2>"
               "<p>This is a minimal framework for 3D mesh manipulation, which can be "
               "extended and used as a base for various projects, for example "
               "prototypes for fairing algorithms, or even displaying/modifying "
               "parametric surfaces, etc.</p>"
               "<p>The following hotkeys are available:</p>"
               "<ul>"
               "<li>&nbsp;M: Toggle mean map</li>"
               "<li>&nbsp;S: Toggle solid (filled polygon) visualization</li>"
               "<li>&nbsp;W: Toggle wireframe visualization</li>"
               "<li>&nbsp;F: Fair mesh</li>"
               "</ul>"
               "<p>There is also a simple selection and movement interface, enabled "
               "only when the wireframe is displayed: a mesh vertex can be selected "
               "by shift-clicking, and it can be moved by shift-dragging. "
               "The vertex always remains on the original mesh.</p>"
               "<p>This is evidently of little practical use; it serves "
               "only to demonstrate the selection and movement process.</p>"
               "<p>Note that libQGLViewer is furnished with a lot of useful features, "
               "such as storing/loading view positions, or saving screenshots. "
               "OpenMesh also has a nice collection of tools for mesh manipulation: "
               "decimation, subdivision, smoothing, etc. These can provide "
               "good comparisons to the methods you implement.</p>"
               "<p>This software can be used as a sample GUI base for handling "
               "parametric or procedural surfaces, as well. The power of "
               "Qt and libQGLViewer makes it easy to set up a prototype application. "
               "Feel free to modify and explore!</p>"
               "<p align=\"right\">Peter Salvi</p>");
  return text;
}
