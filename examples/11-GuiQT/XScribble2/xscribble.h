#ifndef XSCRIBBLE_H
#define XSCRIBBLE_H

#include <QtGui/QMainWindow>
#include "ui_xscribble.h"

class XScribble : public QMainWindow
{
	Q_OBJECT

public:
	XScribble(QWidget *parent = 0, Qt::WFlags flags = 0);
	~XScribble();
private slots:
	void open();
    void openRecentFile();
	void save();	
	void updateStatusBar();
protected:
	void closeEvent(QCloseEvent *event); // overriding ¸²¸Ç
private:
	void createActions();
	void createMenus();
	void createStatusBar();	
	void readSettings();
	void writeSettings();
	void saveSettings();
	bool maybeSave();
	bool saveFile(const QByteArray &fileFormat);
	void updateRecentFileActions();
	Ui::XScribbleClass ui;
	QList<QAction *> m_saveAsActs;
	QMenu* m_saveAsMenu;
	QStringList m_recentFiles;
	enum { MaxRecentFiles = 5 };
	QAction *m_recentFileActions[MaxRecentFiles];
    QAction *m_separatorAction;
	QLabel* m_statusLabel;
};

#endif // XSCRIBBLE_H
