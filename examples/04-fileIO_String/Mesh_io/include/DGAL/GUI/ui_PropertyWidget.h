/********************************************************************************
** Form generated from reading ui file 'PropertyWidget.ui'
**
** Created: Thu Mar 12 09:53:40 2009
**      by: Qt User Interface Compiler version 4.3.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_PROPERTYWIDGET_H
#define UI_PROPERTYWIDGET_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

class Ui_PropertyWidgetClass
{
public:
    QVBoxLayout *vboxLayout;
    QGridLayout *gridLayout;
    QLabel *label_2;
    QComboBox *comboBoxRepresentation;
    QLabel *label;
    QDoubleSpinBox *spinOpacity;
    QSpacerItem *spacerItem;

    void setupUi(QWidget *PropertyWidgetClass)
    {
    if (PropertyWidgetClass->objectName().isEmpty())
        PropertyWidgetClass->setObjectName(QString::fromUtf8("PropertyWidgetClass"));
    QSize size(181, 201);
    size = size.expandedTo(PropertyWidgetClass->minimumSizeHint());
    PropertyWidgetClass->resize(size);
    vboxLayout = new QVBoxLayout(PropertyWidgetClass);
    vboxLayout->setSpacing(6);
    vboxLayout->setMargin(11);
    vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
    gridLayout = new QGridLayout();
    gridLayout->setSpacing(6);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
    label_2 = new QLabel(PropertyWidgetClass);
    label_2->setObjectName(QString::fromUtf8("label_2"));

    gridLayout->addWidget(label_2, 0, 0, 1, 1);

    comboBoxRepresentation = new QComboBox(PropertyWidgetClass);
    comboBoxRepresentation->setObjectName(QString::fromUtf8("comboBoxRepresentation"));

    gridLayout->addWidget(comboBoxRepresentation, 0, 1, 1, 1);

    label = new QLabel(PropertyWidgetClass);
    label->setObjectName(QString::fromUtf8("label"));

    gridLayout->addWidget(label, 1, 0, 1, 1);

    spinOpacity = new QDoubleSpinBox(PropertyWidgetClass);
    spinOpacity->setObjectName(QString::fromUtf8("spinOpacity"));
    spinOpacity->setDecimals(1);
    spinOpacity->setMaximum(1);
    spinOpacity->setSingleStep(0.1);
    spinOpacity->setValue(1);

    gridLayout->addWidget(spinOpacity, 1, 1, 1, 1);


    vboxLayout->addLayout(gridLayout);

    spacerItem = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

    vboxLayout->addItem(spacerItem);


    retranslateUi(PropertyWidgetClass);

    QMetaObject::connectSlotsByName(PropertyWidgetClass);
    } // setupUi

    void retranslateUi(QWidget *PropertyWidgetClass)
    {
    PropertyWidgetClass->setWindowTitle(QApplication::translate("PropertyWidgetClass", "PropertyWidget", 0, QApplication::UnicodeUTF8));
    label_2->setText(QApplication::translate("PropertyWidgetClass", "Representation", 0, QApplication::UnicodeUTF8));
    comboBoxRepresentation->clear();
    comboBoxRepresentation->insertItems(0, QStringList()
     << QApplication::translate("PropertyWidgetClass", "Point", 0, QApplication::UnicodeUTF8)
     << QApplication::translate("PropertyWidgetClass", "Wireframe", 0, QApplication::UnicodeUTF8)
     << QApplication::translate("PropertyWidgetClass", "Surface", 0, QApplication::UnicodeUTF8)
    );
    label->setText(QApplication::translate("PropertyWidgetClass", "Opacity", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(PropertyWidgetClass);
    } // retranslateUi

};

namespace Ui {
    class PropertyWidgetClass: public Ui_PropertyWidgetClass {};
} // namespace Ui

#endif // UI_PROPERTYWIDGET_H
