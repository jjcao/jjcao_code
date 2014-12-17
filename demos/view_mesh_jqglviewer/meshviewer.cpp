#include "meshviewer.h"
#include <QtGui/qfiledialog.h>
#include <QtGui/qstringlistmodel.h>

void MeshViewer::logging(QString& str)
{
	QStringListModel* model = dynamic_cast<QStringListModel*>(ui.listViewPickedVerts->model());
	//int row = ui.listViewPickedVerts->currentIndex().row();
	int row = model->rowCount();
	bool result = model->insertRows(row, 1);
	QModelIndex index = model->index(row);
	result = model->setData (index, str) ;
	ui.listViewPickedVerts->setCurrentIndex(index);
}
void MeshViewer::setupActions()
{
	connect( ui.actionOpen, SIGNAL(triggered()), this, SLOT(open()));
	connect( ui.actionInvert_Face, SIGNAL(triggered()), ui.viewer, SLOT(invertFace()));
	connect( ui.actionInvert_Normal, SIGNAL(triggered()), ui.viewer, SLOT(invertNormal()));
	connect( ui.actionShow_Whole, SIGNAL(triggered()), ui.viewer, SLOT(showWhole()));
	connect( ui.actionShow_Scalar, SIGNAL(triggered()), ui.viewer, SLOT(showScalar()));
	connect( ui.actionShortest_Dist, SIGNAL(triggered()), ui.viewer, SLOT(computeShortestDistance()));

	connect( ui.actionClear_selected, SIGNAL(triggered()), ui.viewer, SLOT(clearSelectedPoints()));
	connect( ui.actionInvert_selected, SIGNAL(triggered()), ui.viewer, SLOT(invertSelectedPoints()));
	connect( ui.actionSave_selected, SIGNAL(triggered()), ui.viewer, SLOT(saveSelectedPoints()));

	connect(ui.viewer, SIGNAL(vertsPicked(QString&)), this, SLOT(logging(QString&)));
}

void MeshViewer::open(QString& filename)
{
	ui.viewer->openMesh(filename);
	setWindowTitle(tr("%1[*] - %2").arg(QFileInfo(filename).fileName())
									.arg(tr("MeshViewer")));
}

void MeshViewer::open()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath());
	if (!filename.isEmpty())		open(filename);
}

MeshViewer::MeshViewer(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	QStringListModel* model = new QStringListModel(ui.listViewPickedVerts);
	ui.listViewPickedVerts->setModel(model);

	setupActions();
}

MeshViewer::~MeshViewer()
{
}