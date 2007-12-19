#include <QString>
#include <KApplication>
#include <KAboutData>
#include <KMessageBox>
#include <KCmdLineArgs>
#include <KLocalizedString>
#include <KDatePicker>
#include <KCalendarSystem>
 
int main (int argc, char *argv[])
{
    KAboutData aboutData("tutorial1",                  // The program name used internally.
                         0,                            // The message catalog name, use program name if null.
                         ki18n("Tutorial 1"),          // A displayable program name string.
                         "1.0",                        // The program version string.
                         ki18n("KMessageBox popup"),   // A short description of what the program does.
                         KAboutData::License_GPL,      // License identifier
                         ki18n("(c) 2007"),            // Copyright Statement
                         ki18n("Some text..."),        // Some free form text, that can contain any kind of information.
                         "http://tutorial.com",        // The program homepage string.
                         "submit@bugs.kde.org");       // The bug report email address string.
 
    KCmdLineArgs::init( argc, argv, &aboutData );
    KApplication app;

    QString cs;

/*
    KDatePicker* picker1 = new KDatePicker();
    picker1->setWindowTitle( "default" );
    picker1->show();

    cs = "hebrew";
    KDatePicker* picker2 = new KDatePicker();
    picker2->setWindowTitle( cs );
    picker2->setCalendar( cs );
    picker2->show();

    cs = "hijri";
    KDatePicker* picker3 = new KDatePicker();
    picker3->setWindowTitle( cs );
    picker3->setCalendar( cs );
    picker3->show();

    cs = "jalali";
    KDatePicker* picker6 = new KDatePicker();
    picker6->setWindowTitle( cs );
    picker6->setCalendar( cs );
    picker6->show();

    cs = "gregorian";
    KDatePicker* picker4 = new KDatePicker();
    picker4->setWindowTitle( cs );
    picker4->setCalendar( cs );
    picker4->show();
*/

    cs = "chinese";
    KDatePicker* picker5 = new KDatePicker();
    picker5->setWindowTitle( cs );
    picker5->setCalendar( cs );
    picker5->show();

    /*
    KGuiItem guiItem( QString( "Hello" ), QString(),
                      QString( "this is a tooltip" ),
                      QString( "this is a whatsthis" ) );
    KMessageBox::questionYesNo( 0, "Hello World", "Hello", guiItem );
    */

    return app.exec();
}

