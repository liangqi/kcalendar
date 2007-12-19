#include <QString>
#include <QDate>
#include <KCalendarSystem>
#include <kcalendarsystemchinese.h>

#include <QtDebug>
 
int main (int argc, char *argv[])
{
  KCalendarSystem *csc = KCalendarSystem::create( "chinese" );
  //KCalendarSystemChinese csc;

  qDebug() << "calendarType: " << csc->calendarType();
  qDebug() << "epoch: " << csc->epoch().toString();
  qDebug() << "earliestValidDate: " << csc->earliestValidDate().toString();
  qDebug() << "latestValidDate: " << csc->latestValidDate().toString();

  qDebug() << "isValid: Chinese 4341-12-30(1645-1-27)" << csc->isValid( 4341, 12, 30 );
  qDebug() << "isValid: Chinese 4342-1-1(1645-1-28)" << csc->isValid( 4342, 1, 1 );
  qDebug() << "isValid: Chinese 9999-12-30(7303-2-5)" << csc->isValid( 9999, 12, 30 );
  qDebug() << "isValid: Chinese 10000-1-1(7303-2-6)" << csc->isValid( 10000, 1, 1 );

  qDebug() << "isValid: 1645-1-27" << csc->isValid( QDate( 1645, 1, 27 ) );
  qDebug() << "isValid: 1645-1-28" << csc->isValid( QDate( 1645, 1, 28 ) );
  qDebug() << "isValid: 7303-2-5" << csc->isValid( QDate( 7303, 2, 5 ) );
  qDebug() << "isValid: 7303-2-6" << csc->isValid( QDate( 7303, 2, 6 ) );

  //qDebug() << "weekStartDay: " << csc->weekStartDay();
  qDebug() << "weekDayOfPray: " << csc->weekDayOfPray();
  qDebug() << "isLunar: " << csc->isLunar();
  qDebug() << "isLunisolar: " << csc->isLunisolar();
  qDebug() << "isSolar: " << csc->isSolar();
  qDebug() << "isProleptic: " << csc->isProleptic();

  bool sign;
  QDate date;

  sign = csc->setDate( date, 4341, 12, 30 );
  qDebug() << "setDate: Chinese 4341-12-30(1645-1-27) " << sign << ", date: " << date.toString();
  sign = csc->setDate( date, 4342, 1, 1 );
  qDebug() << "setDate: Chinese 4342-1-1(1645-1-28) " << sign << ", date: " << date.toString();
  sign = csc->setDate( date, 9999, 12, 30 );
  qDebug() << "setDate: Chinese 9999-12-30(7303-2-5) " << sign << ", date: " << date.toString();
  sign = csc->setDate( date, 10000, 1, 1 );
  qDebug() << "setDate: Chinese 10000-1-1(7303-2-6) " << sign << ", date: " << date.toString();

  sign = csc->setYMD( date, 4341, 12, 30 );
  qDebug() << "setYMD: 4341-12-30 " << sign << ", date: " << date.toString();
  sign = csc->setYMD( date, 4342, 1, 1 );
  qDebug() << "setYMD: 4342-1-1 " << sign << ", date: " << date.toString();
  sign = csc->setYMD( date, 9999, 12, 30 );
  qDebug() << "setYMD: 9999-12-30 " << sign << ", date: " << date.toString();
  sign = csc->setYMD( date, 10000, 1, 1 );
  qDebug() << "setYMD: 10000-1-1 " << sign << ", date: " << date.toString();

  date = QDate( 1645, 1, 27 );
  qDebug() << "year: Chinese 4341-12-30(1645-1-27): " << csc->year( date );
  qDebug() << "month: Chinese 4341-12-30(1645-1-27): " << csc->month( date );
  qDebug() << "day: Chinese 4341-12-30(1645-1-27): " << csc->day( date );

  date = QDate( 1645, 1, 28 );
  qDebug() << "year: Chinese 4342-1-1(1645-1-28): " << csc->year( date );
  qDebug() << "month: Chinese 4342-1-1(1645-1-28): " << csc->month( date );
  qDebug() << "day: Chinese 4342-1-1(1645-1-28): " << csc->day( date );

  date = QDate( 7303, 2, 5 );
  qDebug() << "year: Chinese 9999-12-30(7303-2-5): " << csc->year( date );
  qDebug() << "month: Chinese 9999-12-30(7303-2-5): " << csc->month( date );
  qDebug() << "day: Chinese 9999-12-30(7303-2-5): " << csc->day( date );

  date = QDate( 7303, 2, 6 );
  qDebug() << "year: Chinese 10000-1-1(7303-2-6): " << csc->year( date );
  qDebug() << "month: Chinese 10000-1-1(7303-2-6): " << csc->month( date );
  qDebug() << "day: Chinese 10000-1-1(7303-2-6): " << csc->day( date );

  QDate date1;

/*
  date = QDate( 2007, 12, 19 );
  qDebug() << "Chinese for 2007-12-19: " << csc->year( date ) << "-" << csc->month( date ) << "-" << csc->day( date );

  date1 = csc->addYears( date, 1 );
  qDebug() << "addYears: 2007-12-19, 1: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addYears( date, 2 );
  qDebug() << "addYears: 2007-12-19, 2: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addYears( date, 3 );
  qDebug() << "addYears: 2007-12-19, 3: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addYears( date, 4 );
  qDebug() << "addYears: 2007-12-19, 4: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addYears( date, 5 );
  qDebug() << "addYears: 2007-12-19, 5: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addYears( date, 10 );
  qDebug() << "addYears: 2007-12-19, 10: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addYears( date, 100 );
  qDebug() << "addYears: 2007-12-19, 100: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addYears( date, 1000 );
  qDebug() << "addYears: 2007-12-19, 1000: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );

  date1 = csc->addMonths( date, 1 );
  qDebug() << "addMonths: 2007-12-19, 1: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addMonths( date, 2 );
  qDebug() << "addMonths: 2007-12-19, 2: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );  date1 = csc->addMonths( date, 1 );
  date1 = csc->addMonths( date, 3 );
  qDebug() << "addMonths: 2007-12-19, 3: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );  date1 = csc->addMonths( date, 1 );
  date1 = csc->addMonths( date, 4 );
  qDebug() << "addMonths: 2007-12-19, 4: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );  date1 = csc->addMonths( date, 1 );
  date1 = csc->addMonths( date, 5 );
  qDebug() << "addMonths: 2007-12-19, 5: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );  date1 = csc->addMonths( date, 1 );
  date1 = csc->addMonths( date, 10 );
  qDebug() << "addMonths: 2007-12-19, 10: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );  date1 = csc->addMonths( date, 1 );
  date1 = csc->addMonths( date, 100 );
  qDebug() << "addMonths: 2007-12-19, 100: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );  date1 = csc->addMonths( date, 1 );
  date1 = csc->addMonths( date, 1000 );
  qDebug() << "addMonths: 2007-12-19, 1000: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );  date1 = csc->addMonths( date, 1 );

  date1 = csc->addDays( date, 1 );
  qDebug() << "addDays: 2007-12-19, 1: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addDays( date, 2 );
  qDebug() << "addDays: 2007-12-19, 2: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addDays( date, 3 );
  qDebug() << "addDays: 2007-12-19, 3: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addDays( date, 4 );
  qDebug() << "addDays: 2007-12-19, 4: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addDays( date, 5 );
  qDebug() << "addDays: 2007-12-19, 5: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addDays( date, 10 );
  qDebug() << "addDays: 2007-12-19, 10: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addDays( date, 100 );
  qDebug() << "addDays: 2007-12-19, 100: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
  date1 = csc->addDays( date, 1000 );
  qDebug() << "addDays: 2007-12-19, 1000: " << date1.toString() << ": "  << csc->year( date1 ) << "-" << csc->month( date1 ) << "-" << csc->day( date1 );
*/

  int yearNum;
  date = QDate( 2007, 12, 19 );
  qDebug() << "Chinese for 2007-12-19: " << csc->year( date ) << "-" << csc->month( date ) << "-" << csc->day( date );

  qDebug() << "monthsInYear: " << csc->monthsInYear( date );
  qDebug() << "weeksInYear: " << csc->weeksInYear( date );

  qDebug() << "weeksInYear: " << csc->year( date ) << ": " << csc->weeksInYear( csc->year( date ) );

  qDebug() << "daysInYear: " << csc->daysInYear( date );
  qDebug() << "daysInMonth: " << csc->daysInMonth( date );
  qDebug() << "daysInWeek: " << csc->daysInWeek( date );
  qDebug() << "dayOfYear: " << csc->dayOfYear( date );
  qDebug() << "dayOfWeek: " << csc->dayOfWeek( date );

  qDebug() << "weekNumber: " << csc->weekNumber( date, &yearNum );
  qDebug() << "weekNumber: yearNum = " << yearNum;

  qDebug() << "isLeapYear: " << csc->isLeapYear( date );
  qDebug() << "isLeapYear: " << csc->year( date ) << ": " << csc->isLeapYear( csc->year( date ) );

  date = QDate( 2009, 12, 19 );
  qDebug() << "Chinese for 2009-12-19: " << csc->year( date ) << "-" << csc->month( date ) << "-" << csc->day( date );

  qDebug() << "monthsInYear: " << csc->monthsInYear( date );
  qDebug() << "weeksInYear: " << csc->weeksInYear( date );

  qDebug() << "weeksInYear: " << csc->year( date ) << ": " << csc->weeksInYear( csc->year( date ) );

  qDebug() << "daysInYear: " << csc->daysInYear( date );
  qDebug() << "daysInMonth: " << csc->daysInMonth( date );
  qDebug() << "daysInWeek: " << csc->daysInWeek( date );
  qDebug() << "dayOfYear: " << csc->dayOfYear( date );
  qDebug() << "dayOfWeek: " << csc->dayOfWeek( date );

  qDebug() << "weekNumber: " << csc->weekNumber( date, &yearNum );
  qDebug() << "weekNumber: yearNum = " << yearNum;

  qDebug() << "isLeapYear: " << csc->isLeapYear( date );
  qDebug() << "isLeapYear: " << csc->year( date ) << ": " << csc->isLeapYear( csc->year( date ) );

/*
    virtual QString monthName( int month, int year, MonthNameFormat format = LongName ) const;
    virtual QString monthName( const QDate &date, MonthNameFormat format = LongName ) const;

    virtual QString weekDayName( int weekDay, WeekDayNameFormat format = LongDayName ) const;
    virtual QString weekDayName( const QDate &date, WeekDayNameFormat format = LongDayName ) const;

    virtual QString yearString( const QDate & pDate, StringFormat format = LongFormat ) const;
    virtual QString monthString( const QDate &pDate, StringFormat format = LongFormat ) const;
    virtual QString dayString( const QDate &pDate, StringFormat format = LongFormat ) const;

    virtual int yearStringToInteger( const QString &sNum, int &iLength ) const;
    virtual int monthStringToInteger( const QString &sNum, int &iLength ) const;
    virtual int dayStringToInteger( const QString &sNum, int &iLength ) const;

    virtual QString formatDate( const QDate &date, KLocale::DateFormat format = KLocale::LongDate ) const;

    virtual QDate readDate( const QString &str, bool *ok = 0 ) const;
    virtual QDate readDate( const QString &intstr, const QString &fmt, bool *ok = 0 ) const;
    virtual QDate readDate( const QString &str, KLocale::ReadDateFlags flags, bool *ok = 0 ) const;

    virtual int weekStartDay() const;
    virtual int weekDayOfPray () const;

    virtual bool isLunar() const;
    virtual bool isLunisolar() const;
    virtual bool isSolar() const;
    virtual bool isProleptic() const;
*/

  return 0;
}

