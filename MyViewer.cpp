#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <QtGui/QKeyEvent>

#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MyViewer.h"

#ifdef _WIN32
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_BGRA 0x80E1
#endif

MyViewer::MyViewer(QWidget *parent) :
  QGLViewer(parent), resolution(32), isoline_resolution(10),
  mean_min(0.0), mean_max(0.0), mean_cutoff_ratio(0.05),
  gaussian_min(0.0), gaussian_max(0.0), gaussian_cutoff_ratio(0.05),
  show_control_points(false), show_boundaries(false), show_isolines(false),
  show_solid(true), show_wireframe(false),
  visualization(Visualization::PLAIN), slicing_dir(0, 0, 1), slicing_scaling(1),
  hidden(0), hidden_acc(0), last_filename("")
{
}

MyViewer::~MyViewer() {
  glDeleteTextures(1, &isophote_texture);
  glDeleteTextures(1, &environment_texture);
  glDeleteTextures(1, &slicing_texture);
}

void MyViewer::updateCurvatureMinMax() {
  if (meshes.empty())
    return;

  std::vector<double> mean, gaussian;
  for (const auto &mesh : meshes)
    for (auto v : mesh.vertices()) {
      mean.push_back(mesh.data(v).mean);
      gaussian.push_back(mesh.data(v).gaussian);
    }
  size_t n = mean.size();

  std::sort(mean.begin(), mean.end());
  std::sort(gaussian.begin(), gaussian.end());
  size_t k = (double)n * mean_cutoff_ratio;
  mean_min = std::min(mean[k ? k-1 : 0], 0.0);
  mean_max = std::max(mean[k ? n-k : n-1], 0.0);
  k = (double)n * gaussian_cutoff_ratio;
  gaussian_min = std::min(gaussian[k ? k-1 : 0], 0.0);
  gaussian_max = std::max(gaussian[k ? n-k : n-1], 0.0);
}

static Vec HSV2RGB(Vec hsv) {
  // As in Wikipedia
  double c = hsv[2] * hsv[1];
  double h = hsv[0] / 60;
  double x = c * (1 - std::abs(std::fmod(h, 2) - 1));
  double m = hsv[2] - c;
  Vec rgb(m, m, m);
  if (h <= 1)
    return rgb + Vec(c, x, 0);
  if (h <= 2)
    return rgb + Vec(x, c, 0);
  if (h <= 3)
    return rgb + Vec(0, c, x);
  if (h <= 4)
    return rgb + Vec(0, x, c);
  if (h <= 5)
    return rgb + Vec(x, 0, c);
  if (h <= 6)
    return rgb + Vec(c, 0, x);
  return rgb;
}

static Vec colorMap(double d, double min, double max) {
  double red = 0, green = 120, blue = 240; // Hue
  if (d < 0) {
    double alpha = min ? std::min(d / min, 1.0) : 1.0;
    return HSV2RGB({green * (1 - alpha) + blue * alpha, 1, 1});
  }
  double alpha = max ? std::min(d / max, 1.0) : 1.0;
  return HSV2RGB({green * (1 - alpha) + red * alpha, 1, 1});
}

void MyViewer::updateMesh(bool update_mean_range) {
  meshes.clear();
  for (const auto &s : surfaces)
    meshes.push_back(generateMesh(s));
  if (update_mean_range)
    updateCurvatureMinMax();
}

void MyViewer::setupCamera() {
  // Set camera on the model
  Vector box_min, box_max;
  bool first = true;
  for (const auto &s : surfaces)
    for (const auto &p : s.controlPoints()) {
      Vector v(p[0], p[1], p[2]);
      if (first) {
        first = false;
        box_min = v;
        box_max = v;
      } else {
        box_min.minimize(v);
        box_max.maximize(v);
      }
    }
  camera()->setSceneBoundingBox(Vec(box_min.data()), Vec(box_max.data()));
  camera()->showEntireScene();

  slicing_scaling = 20 / (box_max - box_min).max();

  update();
}

bool MyViewer::openQDS(std::string filename, bool update_view) {
  surfaces.clear();
  try {
    std::ifstream f(filename.c_str());
    f.exceptions(std::ios::failbit | std::ios::badbit);
    size_t n;
    f >> n;
    for (size_t i = 0; i < n; ++i) {
      size_t du, dv, n_ku, n_kv;
      double x, y, z;
      Geometry::DoubleVector ku, kv;
      Geometry::PointVector cpts;
      f >> du >> dv;
      f >> n_ku;
      ku.resize(n_ku);
      for (size_t j = 0; j < n_ku; ++j)
        f >> ku[j];
      f >> n_kv;
      kv.resize(n_kv);
      for (size_t j = 0; j < n_kv; ++j)
        f >> kv[j];
      size_t n_cpts = (n_ku - du - 1) * (n_kv - dv - 1);
      for (size_t j = 0; j < n_cpts; ++j) {
        f >> x >> y >> z;
        cpts.emplace_back(x, y, z);
      }
      surfaces.emplace_back(du, dv, ku, kv, cpts);
    }
  } catch(std::ifstream::failure &) {
    return false;
  }
  last_filename = filename;
  updateMesh(update_view);
  if (update_view)
    setupCamera();
  return true;
}

void MyViewer::init() {
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
  glClearColor(0.33, 0.33, 0.43, 1.0);

  QImage img(":/isophotes.png");
  glGenTextures(1, &isophote_texture);
  glBindTexture(GL_TEXTURE_2D, isophote_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img.width(), img.height(), 0, GL_BGRA,
               GL_UNSIGNED_BYTE, img.convertToFormat(QImage::Format_ARGB32).bits());

  QImage img2(":/environment.png");
  glGenTextures(1, &environment_texture);
  glBindTexture(GL_TEXTURE_2D, environment_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img2.width(), img2.height(), 0, GL_BGRA,
               GL_UNSIGNED_BYTE, img2.convertToFormat(QImage::Format_ARGB32).bits());

  glGenTextures(1, &slicing_texture);
  glBindTexture(GL_TEXTURE_1D, slicing_texture);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  static const unsigned char slicing_img[] = { 0b11111111, 0b00011100 };
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 2, 0, GL_RGB, GL_UNSIGNED_BYTE_3_3_2, &slicing_img);
}

void MyViewer::draw() {
  if (show_control_points)
    for (size_t i = 0; i < surfaces.size(); ++i)
      if (hidden != i + 1)
        drawControlNet(surfaces[i]);

  if (show_boundaries)
    for (size_t i = 0; i < surfaces.size(); ++i)
      if (hidden != i + 1)
      drawBoundaries(surfaces[i]);

  if (show_isolines)
    for (size_t i = 0; i < surfaces.size(); ++i)
      if (hidden != i + 1)
        drawIsolines(surfaces[i]);

  glPolygonMode(GL_FRONT_AND_BACK, !show_solid && show_wireframe ? GL_LINE : GL_FILL);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1, 1);

  if (show_solid || show_wireframe) {
    if (visualization == Visualization::PLAIN)
      glColor3d(1.0, 1.0, 1.0);
    else if (visualization == Visualization::ISOPHOTES) {
      glBindTexture(GL_TEXTURE_2D, current_isophote_texture);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
      glEnable(GL_TEXTURE_2D);
      glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
      glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
      glEnable(GL_TEXTURE_GEN_S);
      glEnable(GL_TEXTURE_GEN_T);
    } else if (visualization == Visualization::SLICING) {
      glBindTexture(GL_TEXTURE_1D, slicing_texture);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
      glEnable(GL_TEXTURE_1D);
    }
    for (size_t i = 0; i < meshes.size(); ++i)
      if (hidden != i + 1) {
        const auto &mesh = meshes[i];
        for (auto f : mesh.faces()) {
          glBegin(GL_POLYGON);
          for (auto v : mesh.fv_range(f)) {
            if (visualization == Visualization::MEAN)
              glColor3dv(colorMap(mesh.data(v).mean, mean_min, mean_max));
            else if (visualization == Visualization::GAUSSIAN)
              glColor3dv(colorMap(mesh.data(v).gaussian, gaussian_min, gaussian_max));
            else if (visualization == Visualization::SLICING)
              glTexCoord1d(mesh.point(v) | slicing_dir * slicing_scaling);
            glNormal3dv(mesh.normal(v).data());
            glVertex3dv(mesh.point(v).data());
          }
          glEnd();
        }
      }
    if (visualization == Visualization::ISOPHOTES) {
      glDisable(GL_TEXTURE_GEN_S);
      glDisable(GL_TEXTURE_GEN_T);
      glDisable(GL_TEXTURE_2D);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    } else if (visualization == Visualization::SLICING) {
      glDisable(GL_TEXTURE_1D);
    }
  }

  if (show_solid && show_wireframe) {
    glPolygonMode(GL_FRONT, GL_LINE);
    glColor3d(0.0, 0.0, 0.0);
    glDisable(GL_LIGHTING);
    for (size_t i = 0; i < meshes.size(); ++i)
      if (hidden != i + 1) {
        const auto &mesh = meshes[i];
        for (auto f : mesh.faces()) {
          glBegin(GL_POLYGON);
          for (auto v : mesh.fv_range(f))
            glVertex3dv(mesh.point(v).data());
          glEnd();
        }
      }
    glEnable(GL_LIGHTING);
  }
}

void MyViewer::drawControlNet(const Geometry::BSSurface &surface) const {
  glDisable(GL_LIGHTING);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glLineWidth(3.0);
  glColor3d(0.3, 0.3, 1.0);
  auto n_cpts = surface.numControlPoints();
  for (size_t k = 0; k < 2; ++k)
    for (size_t i = 0; i < n_cpts[k]; ++i) {
      glBegin(GL_LINE_STRIP);
      for (size_t j = 0; j < n_cpts[1-k]; ++j)
        glVertex3dv(surface.controlPoint(k ? i : j, k ? j : i).data());
      glEnd();
    }
  glLineWidth(1.0);
  glDisable(GL_BLEND);
  glDisable(GL_LINE_SMOOTH);
  glEnable(GL_LIGHTING);
}

void MyViewer::drawBoundaries(const Geometry::BSSurface &surface) const {
  glDisable(GL_LIGHTING);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glLineWidth(5.0);
  glColor3d(0.0, 0.0, 0.0);
  for (size_t k = 0; k < 4; ++k) {
    glBegin(GL_LINE_STRIP);
    for (size_t i = 0; i < resolution; ++i) {
      double u = (double)i / (resolution - 1), v;
      if (k % 2 == 0)
        u = surface.basisU().low() * (1 - u) + surface.basisU().high() * u;
      else
        u = surface.basisV().low() * (1 - u) + surface.basisV().high() * u;
      switch (k) {
      case 0: v = surface.basisV().low(); break;
      case 1: v = surface.basisU().low(); break;
      case 2: v = surface.basisV().high(); break;
      case 3: v = surface.basisU().high(); break;
      default: ;
      }
      auto p = k % 2 == 0 ? surface.eval(u, v) : surface.eval(v, u);
      glVertex3dv(p.data());
    }
    glEnd();
  }
  glLineWidth(1.0);
  glDisable(GL_BLEND);
  glDisable(GL_LINE_SMOOTH);
  glEnable(GL_LIGHTING);
}

void MyViewer::drawIsolines(const Geometry::BSSurface &surface) const {
  glDisable(GL_LIGHTING);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glLineWidth(3.0);
  glColor3d(0.3, 0.7, 0.7);
  for (size_t k = 0; k < 2; ++k) {
    for (size_t i = 0; i <= isoline_resolution; ++i) {
      double u = (double)i / isoline_resolution;
      if (k == 0)
        u = surface.basisU().low() * (1 - u) + surface.basisU().high() * u;
      else
        u = surface.basisV().low() * (1 - u) + surface.basisV().high() * u;
      glBegin(GL_LINE_STRIP);
      for (size_t j = 0; j < resolution; ++j) {
        double v = (double)j / (resolution - 1);
        if (k == 0)
          v = surface.basisV().low() * (1 - v) + surface.basisV().high() * v;
        else
          v = surface.basisU().low() * (1 - v) + surface.basisU().high() * v;
        auto p = k == 0 ? surface.eval(u, v) : surface.eval(v, u);
        glVertex3dv(p.data());
      }
      glEnd();
    }
  }
  glLineWidth(1.0);
  glDisable(GL_BLEND);
  glDisable(GL_LINE_SMOOTH);
  glEnable(GL_LIGHTING);
}

void MyViewer::keyPressEvent(QKeyEvent *e) {
  if (e->key() == Qt::Key_Question) {
    help();
    return;
  }
  if (e->modifiers() == Qt::NoModifier)
    switch (e->key()) {
    case Qt::Key_R:
      openQDS(last_filename, false);
      update();
      break;
    case Qt::Key_O:
      if (camera()->type() == qglviewer::Camera::PERSPECTIVE)
        camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
      else
        camera()->setType(qglviewer::Camera::PERSPECTIVE);
      update();
      break;
    case Qt::Key_P:
      visualization = Visualization::PLAIN;
      update();
      break;
    case Qt::Key_G:
      visualization = Visualization::GAUSSIAN;
      update();
      break;
    case Qt::Key_M:
      visualization = Visualization::MEAN;
      update();
      break;
    case Qt::Key_L:
      visualization = Visualization::SLICING;
      update();
      break;
    case Qt::Key_I:
      visualization = Visualization::ISOPHOTES;
      current_isophote_texture = isophote_texture;
      update();
      break;
    case Qt::Key_E:
      visualization = Visualization::ISOPHOTES;
      current_isophote_texture = environment_texture;
      update();
      break;
    case Qt::Key_C:
      show_control_points = !show_control_points;
      update();
      break;
    case Qt::Key_B:
      show_boundaries = !show_boundaries;
      update();
      break;
    case Qt::Key_N:
      show_isolines = !show_isolines;
      update();
      break;
    case Qt::Key_S:
      show_solid = !show_solid;
      update();
      break;
    case Qt::Key_W:
      show_wireframe = !show_wireframe;
      update();
      break;
    case Qt::Key_0: hidden_acc = hidden_acc * 10 + 0; break;
    case Qt::Key_1: hidden_acc = hidden_acc * 10 + 1; break;
    case Qt::Key_2: hidden_acc = hidden_acc * 10 + 2; break;
    case Qt::Key_3: hidden_acc = hidden_acc * 10 + 3; break;
    case Qt::Key_4: hidden_acc = hidden_acc * 10 + 4; break;
    case Qt::Key_5: hidden_acc = hidden_acc * 10 + 5; break;
    case Qt::Key_6: hidden_acc = hidden_acc * 10 + 6; break;
    case Qt::Key_7: hidden_acc = hidden_acc * 10 + 7; break;
    case Qt::Key_8: hidden_acc = hidden_acc * 10 + 8; break;
    case Qt::Key_9: hidden_acc = hidden_acc * 10 + 9; break;
    case Qt::Key_H:
      if (hidden == 0)
        hidden = hidden_acc;
      else
        hidden = 0;
      hidden_acc = 0;
      update();
      break;
    default:
      QGLViewer::keyPressEvent(e);
    }
  else if (e->modifiers() == Qt::KeypadModifier)
    switch (e->key()) {
    case Qt::Key_Plus:
      if (visualization == Visualization::SLICING)
        slicing_scaling *= 2;
      else if (show_isolines)
        isoline_resolution *= 2;
      else {
        resolution *= 2;
        updateMesh();
      }
      update();
      break;
    case Qt::Key_Minus:
      if (visualization == Visualization::SLICING)
        slicing_scaling /= 2;
      else if (show_isolines) {
        if (isoline_resolution > 5)
          isoline_resolution /= 2;
      } else if (resolution > 4) {
        resolution /= 2;
        updateMesh();
      }
      update();
      break;
    case Qt::Key_Asterisk:
      slicing_dir = Vector(static_cast<double *>(camera()->viewDirection()));
      update();
      break;
    } else
    QGLViewer::keyPressEvent(e);
}

MyViewer::MyMesh MyViewer::generateMesh(const Geometry::BSSurface &surface) {
  MyMesh mesh;
  mesh.request_vertex_normals();

  std::vector<MyMesh::VertexHandle> handles, tri;
  std::vector<std::pair<MyMesh::VertexHandle, Vector>> normals;

  Geometry::VectorMatrix der;
  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / (double)(resolution - 1);
    u = surface.basisU().low() * (1 - u) + surface.basisU().high() * u;
    for (size_t j = 0; j < resolution; ++j) {
      double v = (double)j / (double)(resolution - 1);
      v = surface.basisV().low() * (1 - v) + surface.basisV().high() * v;
      auto p = surface.eval(u, v, 2, der);
      auto h = mesh.add_vertex(Vector(p.data()));
      auto &Su = der[1][0];
      auto &Sv = der[0][1];
      auto &Suu = der[2][0];
      auto &Suv = der[1][1];
      auto &Svv = der[0][2];
      auto n = (Su ^ Sv).normalize();
      auto E = Su.normSqr();
      auto F = Su * Sv;
      auto G = Sv.normSqr();
      auto L = Suu * n;
      auto M = Suv * n;
      auto N = Svv * n;
      mesh.set_normal(h, Vector((-n).data()));
      mesh.data(h).mean = (N * E - 2 * M * F + L * G) / (2 * (E * G - F * F));
      mesh.data(h).gaussian = (L * N - M * M) / (E * G - F * F);
      handles.push_back(h);
    }
  }
  for (size_t i = 0; i < resolution - 1; ++i)
    for (size_t j = 0; j < resolution - 1; ++j) {
      tri.clear();
      tri.push_back(handles[i * resolution + j]);
      tri.push_back(handles[i * resolution + j + 1]);
      tri.push_back(handles[(i + 1) * resolution + j]);
      mesh.add_face(tri);
      tri.clear();
      tri.push_back(handles[(i + 1) * resolution + j]);
      tri.push_back(handles[i * resolution + j + 1]);
      tri.push_back(handles[(i + 1) * resolution + j + 1]);
      mesh.add_face(tri);
    }
  return mesh;
}

QString MyViewer::helpString() const {
  QString text("<h2>QDS Viewer</h2>"
               "<p>This is a quad patch viewer.</p>"
               "<p>The following hotkeys are available:</p>"
               "<ul>"
               "<li>&nbsp;R: Reload model</li>"
               "<li>&nbsp;O: Toggle orthographic projection</li>"
               "<li>&nbsp;P: Set plain map (no coloring)</li>"
               "<li>&nbsp;G: Set Gaussian curvature map</li>"
               "<li>&nbsp;M: Set mean curvature map</li>"
               "<li>&nbsp;L: Set slicing map<ul>"
               "<li>&nbsp;+: Increase resolution or isoline/slicing density</li>"
               "<li>&nbsp;-: Decrease resolution or isoline/slicing density</li>"
               "<li>&nbsp;*: Set slicing direction to view</li></ul></li>"
               "<li>&nbsp;I: Set isophote line map</li>"
               "<li>&nbsp;E: Set environment texture</li>"
               "<li>&nbsp;C: Toggle control polygon visualization</li>"
               "<li>&nbsp;B: Toggle boundary curve visualization</li>"
               "<li>&nbsp;N: Toggle isoline curve visualization</li>"
               "<li>&nbsp;S: Toggle solid (filled polygon) visualization</li>"
               "<li>&nbsp;W: Toggle wireframe visualization</li>"
               "<li>&nbsp;&lt;n&gt;H: Hide surface #n (H shows it again)</li>"
               "<li>&nbsp;?: This help</li>"
               "</ul>"
               "<p align=\"right\">Peter Salvi</p>");
  return text;
}
