// -*- mode: c++ -*-
#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <geometry.hh>

using qglviewer::Vec;

class MyViewer : public QGLViewer {
  Q_OBJECT

public:
  explicit MyViewer(QWidget *parent);
  virtual ~MyViewer();

  inline double getCutoffRatio() const;
  inline void setCutoffRatio(double ratio);
  inline double getMeanMin() const;
  inline void setMeanMin(double min);
  inline double getMeanMax() const;
  inline void setMeanMax(double max);
  inline const double *getSlicingDir() const;
  inline void setSlicingDir(double x, double y, double z);
  inline double getSlicingScaling() const;
  inline void setSlicingScaling(double scaling);
  bool openQDS(std::string filename, bool update_view = true);

signals:
  void startComputation(QString message);
  void midComputation(int percent);
  void endComputation();

protected:
  virtual void init() override;
  virtual void draw() override;
  virtual void keyPressEvent(QKeyEvent *e) override;
  virtual QString helpString() const override;

private:
  struct MyTraits : public OpenMesh::DefaultTraits {
    using Point  = OpenMesh::Vec3d; // the default would be Vec3f
    using Normal = OpenMesh::Vec3d;
    VertexTraits {
      double mean;              // approximated mean curvature
    };
  };
  using MyMesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;
  using Vector = OpenMesh::VectorT<double,3>;

  // Mesh
  void updateMesh(bool update_mean_range = true);
  void updateMeanMinMax();

  MyMesh generateMesh(const Geometry::BSSurface &surface);

  // Visualization
  void setupCamera();
  Vec meanMapColor(double d) const;
  void drawControlNet(const Geometry::BSSurface &surface) const;
  void drawBoundaries(const Geometry::BSSurface &surface) const;
  void drawIsolines(const Geometry::BSSurface &surface) const;

  //////////////////////
  // Member variables //
  //////////////////////

  std::vector<Geometry::BSSurface> surfaces;
  std::vector<MyMesh> meshes;
  size_t resolution, isoline_resolution;

  // Visualization
  double mean_min, mean_max, cutoff_ratio;
  bool show_control_points, show_boundaries, show_isolines, show_solid, show_wireframe;
  enum class Visualization { PLAIN, MEAN, SLICING, ISOPHOTES } visualization;
  GLuint isophote_texture, environment_texture, current_isophote_texture, slicing_texture;
  Vector slicing_dir;
  double slicing_scaling;
  std::string last_filename;
};

#include "MyViewer.hpp"
