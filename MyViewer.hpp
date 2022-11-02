double MyViewer::getMeanCutoffRatio() const {
  return mean_cutoff_ratio;
}

void MyViewer::setMeanCutoffRatio(double ratio) {
  mean_cutoff_ratio = ratio;
  updateCurvatureMinMax();
}

double MyViewer::getGaussianCutoffRatio() const {
  return gaussian_cutoff_ratio;
}

void MyViewer::setGaussianCutoffRatio(double ratio) {
  gaussian_cutoff_ratio = ratio;
  updateCurvatureMinMax();
}

double MyViewer::getMeanMin() const {
  return mean_min;
}

void MyViewer::setMeanMin(double min) {
  mean_min = min;
}

double MyViewer::getMeanMax() const {
  return mean_max;
}

void MyViewer::setMeanMax(double max) {
  mean_max = max;
}

double MyViewer::getGaussianMin() const {
  return gaussian_min;
}

void MyViewer::setGaussianMin(double min) {
  gaussian_min = min;
}

double MyViewer::getGaussianMax() const {
  return gaussian_max;
}

void MyViewer::setGaussianMax(double max) {
  gaussian_max = max;
}

const double *MyViewer::getSlicingDir() const {
  return slicing_dir.data();
}

void MyViewer::setSlicingDir(double x, double y, double z) {
  slicing_dir = Vector(x, y, z).normalized();
}

double MyViewer::getSlicingScaling() const {
  return slicing_scaling;
}

void MyViewer::setSlicingScaling(double scaling) {
  slicing_scaling = scaling;
}

