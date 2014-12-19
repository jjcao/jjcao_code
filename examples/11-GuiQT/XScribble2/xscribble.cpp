#include "xscribble.h"
#include <QtGui>
#include <algorithm>

void XScribble::createActions()
{
	connect( ui.actionOpen, SIGNAL(triggered()), this, SLOT(open()));
	foreach (QByteArray format, QImageWriter::supportedImageFormats()) {
		QString text = tr("%1...").arg(QString(format).toUpper());

		QAction *action = new QAction(text, this);
		action->setData(format);
		connect(action, SIGNAL(triggered()), this, SLOT(save()));
		m_saveAsActs.append(action);
	}

	for (int i = 0; i < MaxRecentFiles; ++i) {
        m_recentFileActions[i] = new QAction(this);
        m_recentFileActions[i]->setVisible(false);
        connect(m_recentFileActions[i], SIGNAL(triggered()),
                this, SLOT(openRecentFile()));
    }
	
	ui.actionFreeHand->setData("Free Hand");
	ui.actionLine->setData("Line");
	ui.actionCircle->setData("Circle");
	connect( ui.actionFreeHand, SIGNAL(triggered()), ui.centralWidget, SLOT(flagFree()) );
	connect( ui.actionLine, SIGNAL(triggered()), ui.centralWidget, SLOT(flagLine()) );
	connect( ui.actionCircle, SIGNAL(triggered()), ui.centralWidget, SLOT(flagCircle()) );
	connect( ui.actionSetPenColor, SIGNAL(triggered()), ui.centralWidget, SLOT(setPenColor()) );
}
void XScribble::createStatusBar()
{
    m_statusLabel = new QLabel("Ready!");
    m_statusLabel->setAlignment(Qt::AlignHCenter);
    m_statusLabel->setMinimumSize(m_statusLabel->sizeHint());

	ui.statusBar->addWidget(m_statusLabel);

	connect( ui.actionFreeHand, SIGNAL(triggered()), this, SLOT(updateStatusBar()) );
	connect( ui.actionLine, SIGNAL(triggered()), this, SLOT(updateStatusBar()) );
	connect( ui.actionCircle, SIGNAL(triggered()), this, SLOT(updateStatusBar()) );

    updateStatusBar();
}
void XScribble::updateStatusBar()
{
	 QAction *action = qobject_cast<QAction *>(sender());
    if (action)
	{
		QString tmp = "scribble tools changed to: " + action->data().toString();
		m_statusLabel->setText(tmp);
	}    
}
void XScribble::createMenus()
{
	m_saveAsMenu = new QMenu(tr("&Save As"), this);
    foreach (QAction *action, m_saveAsActs)
        m_saveAsMenu->addAction(action);

	ui.menuFile->insertMenu(ui.actionExit, m_saveAsMenu);

	m_separatorAction = ui.menuFile->addSeparator();
    for (int i = 0; i < MaxRecentFiles; ++i)
        ui.menuFile->addAction(m_recentFileActions[i]);
}
void XScribble::readSettings()
{
    QSettings settings("Ability Co.Ltd.", "XScribble");
    QRect rect = settings.value("geometry",  QRect(200, 200, 400, 400)).toRect();
    move(rect.topLeft());
    resize(rect.size());

    m_recentFiles = settings.value("recentFiles").toStringList();
    updateRecentFileActions();
}
void XScribble::writeSettings()
{
    QSettings settings("Ability Co.Ltd.", "XScribble");

    settings.setValue("geometry", geometry());
    settings.setValue("recentFiles", m_recentFiles);
}
void XScribble::closeEvent(QCloseEvent *event)
{
    if (maybeSave()) {
        writeSettings();
        event->accept();
    } else {
        event->ignore();
    }
}
XScribble::~XScribble()
{
}
XScribble::XScribble(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	createActions();
	createMenus();	
	createStatusBar();

	readSettings();
}
bool XScribble::maybeSave()
{
    if (ui.centralWidget->isModified()) {
       QMessageBox::StandardButton ret;
       ret = QMessageBox::warning(this, tr("Scribble"),
                          tr("The image has been modified.\n"
                             "Do you want to save your changes?"),
                          QMessageBox::Save | QMessageBox::Discard
						| QMessageBox::Cancel);
        if (ret == QMessageBox::Save) {
            return saveFile("png");
        } else if (ret == QMessageBox::Cancel) {
            return false;
        }
    }
    return true;
}
void XScribble::open()
{
	if (maybeSave()) {
        QString fileName = QFileDialog::getOpenFileName(this,
                                   tr("Open File"), QDir::currentPath());
        if (!fileName.isEmpty())
		{
            ui.centralWidget->openImage(fileName);

			m_recentFiles.removeAll(fileName);
			m_recentFiles.prepend(fileName);
			updateRecentFileActions();
			setWindowTitle(tr("%1[*] - %2").arg(QFileInfo(fileName).fileName())
									   .arg(tr("Scribble")));
		}
    }    
}
void XScribble::openRecentFile()
{
    if (maybeSave()) {
        QAction *action = qobject_cast<QAction *>(sender());
        if (action)
		{
			QString fileName(action->data().toString());
			ui.centralWidget->openImage(fileName);

			m_recentFiles.removeAll(fileName);
			m_recentFiles.prepend(fileName);
			updateRecentFileActions();
			setWindowTitle(tr("%1[*] - %2").arg(QFileInfo(fileName).fileName())
									   .arg(tr("Scribble")));
		}
    }
}
void XScribble::save()
{
    QAction *action = qobject_cast<QAction *>(sender());
    QByteArray fileFormat = action->data().toByteArray();
    saveFile(fileFormat);
}
bool XScribble::saveFile(const QByteArray &fileFormat)
{
    QString initialPath = QDir::currentPath() + "/untitled." + fileFormat;

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save As"),
                               initialPath,
                               tr("%1 Files (*.%2);;All Files (*)")
                               .arg(QString(fileFormat.toUpper()))
                               .arg(QString(fileFormat)));
    if (fileName.isEmpty()) {
        return false;
    } else {
        return ui.centralWidget->saveImage(fileName, fileFormat);
    }
}
bool isNotExists(const QString& elem)//functor
{
	bool result = QFile::exists(elem);
	return !result;
}

void XScribble::updateRecentFileActions()
{    
	m_recentFiles.erase( std::remove_if(m_recentFiles.begin(), m_recentFiles.end(), isNotExists), m_recentFiles.end());

	for (int j = 0; j < MaxRecentFiles; ++j) {
        if (j < m_recentFiles.count()) {
            QString text = tr("&%1 %2")
                           .arg(j + 1)
                           .arg(QFileInfo( m_recentFiles[j] ).fileName());
            m_recentFileActions[j]->setText(text);
            m_recentFileActions[j]->setData(m_recentFiles[j]);
            m_recentFileActions[j]->setVisible(true);
        } else {
            m_recentFileActions[j]->setVisible(false);
        }
    }
    m_separatorAction->setVisible(!m_recentFiles.isEmpty());
}