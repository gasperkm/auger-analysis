#include "primary_type.h"

PrimPart::PrimPart()
{
   names = new vector<string>;
   shortnames = new vector<string>;
   masses = new vector<double>;

   names->push_back("Proton");
   names->push_back("Helium");
   names->push_back("Lithium");
   names->push_back("Beryllium");
   names->push_back("Boron");
   names->push_back("Carbon");
   names->push_back("Nitrogen");
   names->push_back("Oxygen");
   names->push_back("Fluorine");
   names->push_back("Neon");
   names->push_back("Natrium");
   names->push_back("Magnesium");
   names->push_back("Aluminium");
   names->push_back("Silicon");
   names->push_back("Phosphorus");
   names->push_back("Sulfur");
   names->push_back("Chlorine");
   names->push_back("Argon");
   names->push_back("Potassium");
   names->push_back("Calcium");
   names->push_back("Scandium");
   names->push_back("Titanium");
   names->push_back("Vanadium");
   names->push_back("Chromium");
   names->push_back("Manganese");
   names->push_back("Iron");
   names->push_back("Other");

   shortnames->push_back("p");
   shortnames->push_back("He");
   shortnames->push_back("Li");
   shortnames->push_back("Be");
   shortnames->push_back("B");
   shortnames->push_back("C");
   shortnames->push_back("N");
   shortnames->push_back("O");
   shortnames->push_back("F");
   shortnames->push_back("Ne");
   shortnames->push_back("Na");
   shortnames->push_back("Mg");
   shortnames->push_back("Al");
   shortnames->push_back("Si");
   shortnames->push_back("P");
   shortnames->push_back("S");
   shortnames->push_back("Cl");
   shortnames->push_back("Ar");
   shortnames->push_back("K");
   shortnames->push_back("Ca");
   shortnames->push_back("Sc");
   shortnames->push_back("Ti");
   shortnames->push_back("V");
   shortnames->push_back("Cr");
   shortnames->push_back("Mn");
   shortnames->push_back("Fe");
   shortnames->push_back("Other");

   masses->push_back(1.008);
   masses->push_back(4.0026);
   masses->push_back(6.94);
   masses->push_back(9.0122);
   masses->push_back(10.81);
   masses->push_back(12.011);
   masses->push_back(14.007);
   masses->push_back(15.999);
   masses->push_back(18.998);
   masses->push_back(20.180);
   masses->push_back(22.990);
   masses->push_back(24.305);
   masses->push_back(26.982);
   masses->push_back(28.085);
   masses->push_back(30.974);
   masses->push_back(32.06);
   masses->push_back(35.45);
   masses->push_back(39.948);
   masses->push_back(39.098);
   masses->push_back(40.078);
   masses->push_back(44.956);
   masses->push_back(47.867);
   masses->push_back(50.942);
   masses->push_back(51.996);
   masses->push_back(54.938);
   masses->push_back(55.845);
   masses->push_back(-1);

   nrelem = names->size();
}

PrimPart::~PrimPart()
{
   delete names;
   delete shortnames;
   delete masses;
}

string PrimPart::GetName(string *inname)
{
   for(int i = 0; i < nrelem; i++)
   {
      if(inname->compare(shortnames->at(i)) == 0)
      {
         stemp = names->at(i);
         return stemp;
      }
   }

   return string("empty");
}

string PrimPart::GetName(int inz)
{
   stemp = names->at(inz);
   return stemp;
}

string PrimPart::GetShortName(string *inname)
{
   stemp = "";

   for(int i = 0; i < nrelem; i++)
   {
      if(inname->compare(names->at(i)) == 0)
      {
         stemp = shortnames->at(i);
         return stemp;
      }
   }

   return string("empty");
}

string PrimPart::GetShortName(int inz)
{
   stemp = shortnames->at(inz);
   return stemp;
}

int PrimPart::GetZ(string *inname)
{
   for(int i = 0; i < nrelem; i++)
   {
      if(inname->compare(names->at(i)) == 0)
      {
         return i;
      }
   }

   for(int i = 0; i < nrelem; i++)
   {
      if(inname->compare(shortnames->at(i)) == 0)
      {
         return i;
      }
   }

   return -1;
}

double PrimPart::GetA(string *inname)
{
   for(int i = 0; i < nrelem; i++)
   {
      if(inname->compare(names->at(i)) == 0)
      {
         dtemp = masses->at(i);
         return dtemp;
      }
   }

   for(int i = 0; i < nrelem; i++)
   {
      if(inname->compare(shortnames->at(i)) == 0)
      {
         dtemp = masses->at(i);
         return dtemp;
      }
   }

   return -1;
}

double PrimPart::GetA(int inz)
{
   dtemp = masses->at(inz);
   return dtemp;
}

int PrimPart::Nr()
{
   return nrelem;
}

/*int PrimPart::Z(string name)
{
   for(int i = 0; i < names->size(); i++)
   {
      if(name.compare(names->at(i)) == 0)
         return i;
   }

   for(int i = 0; i < shortnames->size(); i++)
   {
      if(name.compare(shortnames->at(i)) == 0)
         return i;
   }
}

string PrimPart::Name(string shortname)
{
   for(int i = 0; i < shortnames->size(); i++)
   {
      if(shortname.compare(names->at(i)) == 0)
         return names->at(i);
   }
}*/

/*string PrimPart::ShortName(string name)
{
   cout << "Bu!" << endl;
   cout << "size of shortnames = " << shortnames.size() << endl;
   cout << "shortnames[6] = " << shortnames.at(6) << endl;
   cout << "name = " << name << endl;*/
/*   for(int i = 0; i < shortnames->size(); i++)
   {
      cout << "shortname = " << shortnames->at(i) << endl;
      if(name.compare(shortnames->at(i)) == 0)
         return shortnames->at(i);
   }
}*/

/*string PrimPart::Name(int z)
{
   return names->at(z-1);
}

string PrimPart::ShortName(int z)
{
   return shortnames->at(z-1);
}

double PrimPart::A(string name)
{
   for(int i = 0; i < names->size(); i++)
   {
      if(name.compare(names->at(i)) == 0)
         return masses->at(i);
   }

   for(int i = 0; i < shortnames->size(); i++)
   {
      if(name.compare(shortnames->at(i)) == 0)
         return masses->at(i);
   }
}

double PrimPart::A(int z)
{
   return masses->at(z-1);
}

int PrimPart::N()
{
   return names.size();
}*/
