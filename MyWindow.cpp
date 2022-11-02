#include <memory>

#include <QtWidgets>

#include "MyWindow.h"

MyWindow::MyWindow(QApplication *parent) :
  QMainWindow(), parent(parent), last_directory(".")
{
  setWindowTitle(tr("QDS Viewer"));
  setStatusBar(new QStatusBar);
  progress = new QProgressBar;
  progress->setMinimum(0); progress->setMaximum(100);
  progress->hide();
  statusBar()->addPermanentWidget(progress);

  viewer = new MyViewer(this);
  connect(viewer, SIGNAL(startComputation(QString)), this, SLOT(startComputation(QString)));
  connect(viewer, SIGNAL(midComputation(int)), this, SLOT(midComputation(int)));
  connect(viewer, SIGNAL(endComputation()), this, SLOT(endComputation()));
  setCentralWidget(viewer);

  /////////////////////////
  // Setup actions/menus //
  /////////////////////////

  auto openAction = new QAction(tr("&Open"), this);
  openAction->setShortcut(tr("Ctrl+O"));
  openAction->setStatusTip(tr("Load a model from a file"));
  connect(openAction, SIGNAL(triggered()), this, SLOT(open()));

  auto quitAction = new QAction(tr("&Quit"), this);
  quitAction->setShortcut(tr("Ctrl+Q"));
  quitAction->setStatusTip(tr("Quit the program"));
  connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));

  auto cutoffAction = new QAction(tr("Set &cutoff ratio"), this);
  cutoffAction->setStatusTip(tr("Set curvature map cutoff ratio"));
  connect(cutoffAction, SIGNAL(triggered()), this, SLOT(setCutoff()));

  auto rangeAction = new QAction(tr("Set &range"), this);
  rangeAction->setStatusTip(tr("Set curvature map range"));
  connect(rangeAction, SIGNAL(triggered()), this, SLOT(setRange()));

  auto slicingAction = new QAction(tr("Set &slicing parameters"), this);
  slicingAction->setStatusTip(tr("Set contouring direction and scaling"));
  connect(slicingAction, SIGNAL(triggered()), this, SLOT(setSlicing()));

  auto fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(openAction);
  fileMenu->addAction(quitAction);

  auto visMenu = menuBar()->addMenu(tr("&Visualization"));
  visMenu->addAction(cutoffAction);
  visMenu->addAction(rangeAction);
  visMenu->addAction(slicingAction);
}

MyWindow::~MyWindow() {
}

void MyWindow::open() {
  auto filename =
    QFileDialog::getOpenFileName(this, tr("Open File"), last_directory,
                                 tr("Quad patches (*.qds)"));
  if(filename.isEmpty())
    return;
  last_directory = QFileInfo(filename).absolutePath();

  bool ok = viewer->openQDS(filename.toUtf8().data());

  if (!ok)
    QMessageBox::warning(this, tr("Cannot open file"),
                         tr("Could not open file: ") + filename + ".");
}

void MyWindow::setCutoff() {
  // Memory management options for the dialog:
  // - on the stack (deleted at the end of the function)
  // - on the heap with manual delete or std::unique_ptr 
  // There is also a Qt-controlled automatic deletion by calling
  //     dlg->setAttribute(Qt::WA_DeleteOnClose);
  // ... but this could delete sub-widgets too early.

  auto dlg = std::make_unique<QDialog>(this);
  auto *hb1    = new QHBoxLayout,
       *hb2    = new QHBoxLayout,
       *hb3    = new QHBoxLayout;
  auto *vb     = new QVBoxLayout;
  auto *text1  = new QLabel(tr("Mean cutoff ratio:"));
  auto *sb1    = new QDoubleSpinBox;
  auto *text2  = new QLabel(tr("Gaussian cutoff ratio:"));
  auto *sb2    = new QDoubleSpinBox;
  auto *cancel = new QPushButton(tr("Cancel"));
  auto *ok     = new QPushButton(tr("Ok"));

  sb1->setDecimals(3);
  sb1->setRange(0.001, 0.5);
  sb1->setSingleStep(0.01);
  sb1->setValue(viewer->getMeanCutoffRatio());
  sb2->setDecimals(3);
  sb2->setRange(0.001, 0.5);
  sb2->setSingleStep(0.01);
  sb2->setValue(viewer->getGaussianCutoffRatio());
  connect(cancel, SIGNAL(pressed()), dlg.get(), SLOT(reject()));
  connect(ok,     SIGNAL(pressed()), dlg.get(), SLOT(accept()));
  ok->setDefault(true);

  hb1->addWidget(text1);
  hb1->addWidget(sb1);
  hb2->addWidget(text2);
  hb2->addWidget(sb2);
  hb3->addWidget(cancel);
  hb3->addWidget(ok);
  vb->addLayout(hb1);
  vb->addLayout(hb2);
  vb->addLayout(hb3);

  dlg->setWindowTitle(tr("Set ratio"));
  dlg->setLayout(vb);

  if(dlg->exec() == QDialog::Accepted) {
    viewer->setMeanCutoffRatio(sb1->value());
    viewer->setGaussianCutoffRatio(sb2->value());
    viewer->update();
  }
}

void MyWindow::setRange() {
  QDialog dlg(this);
  auto *grid   = new QGridLayout;
  auto *text1  = new QLabel(tr("Mean min:")),
       *text2  = new QLabel(tr("Mean max:"));
  auto *sb1    = new QDoubleSpinBox,
       *sb2    = new QDoubleSpinBox;
  auto *text3  = new QLabel(tr("Gaussian min:")),
       *text4  = new QLabel(tr("Gaussian max:"));
  auto *sb3    = new QDoubleSpinBox,
       *sb4    = new QDoubleSpinBox;
  auto *cancel = new QPushButton(tr("Cancel"));
  auto *ok     = new QPushButton(tr("Ok"));

  // The range of the spinbox controls the number of displayable digits,
  // so setting it to a large value results in a very wide window.
  double max = 1000.0; // std::numeric_limits<double>::max();
  sb1->setDecimals(5);                 sb2->setDecimals(5);
  sb1->setRange(-max, 0.0);            sb2->setRange(0.0, max);
  sb1->setSingleStep(0.01);            sb2->setSingleStep(0.01);
  sb1->setValue(viewer->getMeanMin()); sb2->setValue(viewer->getMeanMax());
  sb3->setDecimals(5);                 sb4->setDecimals(5);
  sb3->setRange(-max, 0.0);            sb4->setRange(0.0, max);
  sb3->setSingleStep(0.01);            sb4->setSingleStep(0.01);
  sb3->setValue(viewer->getGaussianMin()); sb4->setValue(viewer->getGaussianMax());
  connect(cancel, SIGNAL(pressed()), &dlg, SLOT(reject()));
  connect(ok,     SIGNAL(pressed()), &dlg, SLOT(accept()));
  ok->setDefault(true);

  grid->addWidget( text1, 1, 1, Qt::AlignRight);
  grid->addWidget(   sb1, 1, 2);
  grid->addWidget( text2, 2, 1, Qt::AlignRight);
  grid->addWidget(   sb2, 2, 2);
  grid->addWidget( text3, 3, 1, Qt::AlignRight);
  grid->addWidget(   sb3, 3, 2);
  grid->addWidget( text4, 4, 1, Qt::AlignRight);
  grid->addWidget(   sb4, 4, 2);
  grid->addWidget(cancel, 5, 1);
  grid->addWidget(    ok, 5, 2);

  dlg.setWindowTitle(tr("Set range"));
  dlg.setLayout(grid);

  if(dlg.exec() == QDialog::Accepted) {
    viewer->setMeanMin(sb1->value());
    viewer->setMeanMax(sb2->value());
    viewer->setGaussianMin(sb3->value());
    viewer->setGaussianMax(sb4->value());
    viewer->update();
  }
}

void MyWindow::setSlicing() {
  auto dlg = std::make_unique<QDialog>(this);
  auto *hb1    = new QHBoxLayout,
       *hb2    = new QHBoxLayout,
       *hb3    = new QHBoxLayout;
  auto *vb     = new QVBoxLayout;
  auto *text_v = new QLabel(tr("Slicing direction:"));
  QDoubleSpinBox *sb_v[3];
  auto *text_s = new QLabel(tr("Slicing scaling:"));
  auto *sb_s   = new QDoubleSpinBox;
  auto *cancel = new QPushButton(tr("Cancel"));
  auto *ok     = new QPushButton(tr("Ok"));

  for (int i = 0; i < 3; ++i) {
    sb_v[i] = new QDoubleSpinBox;
    sb_v[i]->setDecimals(3);
    sb_v[i]->setRange(-1, 1);
    sb_v[i]->setSingleStep(0.01);
    sb_v[i]->setValue(viewer->getSlicingDir()[i]);
  }
  sb_s->setDecimals(6);
  sb_s->setRange(0, 10000);
  sb_s->setSingleStep(1);
  sb_s->setValue(viewer->getSlicingScaling());
  connect(cancel, SIGNAL(pressed()), dlg.get(), SLOT(reject()));
  connect(ok,     SIGNAL(pressed()), dlg.get(), SLOT(accept()));
  ok->setDefault(true);

  hb1->addWidget(text_v);
  hb1->addWidget(sb_v[0]); hb1->addWidget(sb_v[1]); hb1->addWidget(sb_v[2]);
  hb2->addWidget(text_s);
  hb2->addWidget(sb_s);
  hb3->addWidget(cancel);
  hb3->addWidget(ok);
  vb->addLayout(hb1);
  vb->addLayout(hb2);
  vb->addLayout(hb3);

  dlg->setWindowTitle(tr("Set slicing"));
  dlg->setLayout(vb);

  if(dlg->exec() == QDialog::Accepted) {
    viewer->setSlicingDir(sb_v[0]->value(), sb_v[1]->value(), sb_v[2]->value());
    viewer->setSlicingScaling(sb_s->value());
    viewer->update();
  }
}

void MyWindow::startComputation(QString message) {
  statusBar()->showMessage(message);
  progress->setValue(0);
  progress->show();
  parent->processEvents(QEventLoop::ExcludeUserInputEvents);
}

void MyWindow::midComputation(int percent) {
  progress->setValue(percent);
  parent->processEvents(QEventLoop::ExcludeUserInputEvents);
}

void MyWindow::endComputation() {
  progress->hide();
  statusBar()->clearMessage();
}
