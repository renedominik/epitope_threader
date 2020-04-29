/////////////////////////////////////////////////////////////////////////
//
//  This file is part of the EpitopeThreader program
//
//  Copyright (C) 2013 by Rene Staritzbichler
//  renedominik@yahoo.de
//
//
//
//
// date:   Dez 16, 2013
/////////////////////////////////////////////////////////////////////////



#ifndef SIDECHAIN_H_
#define SIDECHAIN_H_


class RotamerLibrary
{
private:
	std::map< std::string, std::vector< std::vector< math::Vector3N> > >  m_Rotamers;    //<  Sidechain atom positions

//	std::vector< std::vector< math::Vector3N *> >                         m_MolMirror;   //<  Pointer for each atom in mols to first rotamer


public:
   void
   RandomRotamer( Residue &RES)
   {
	   math::Vector3N
		   pos = RES.Atoms()[1],
		   x = RES.Atoms()[0] - RES.Atoms()[1],
		   y = RES.Atoms()[2] - RES.Atoms()[1],
		   z = CrossProduct( x, y);

	   std::vector< math::Vector3N>
		   *rotamer = &RandomElement( m_Rotamers[ RES.Name()]);

	   // kick out if possible
	   if( rotamer->size() != RES.Atoms().size() - 3)
	   {
		   BuildSidechain( RES);
	   }

	   std::vector< Atom>::iterator res = RES.Atoms().begin() + 3;
	   for( std::vector< math::Vector3N>::iterator rot = rotamer->begin(); rot != rotamer->end(); ++rot, ++res)
	   {
		   res->Position() = pos + (*rot)[0] * x + (*rot)[1] * y + (*rot)[2] * z;
	   }
   }

   void BuildSidechain( Residue &RES)
   {
	   std::ifstream in;
	   Open( __PATH__ + "../database/" + RES.Name() + "_default.txt", in);

	   Close( in);
   }


   void
   ExtractRotamer( const Residue &RES)
   {
	   std::vector< math::Vector3N>
		   v, coo(4);

	   coo[0] = RES.Atoms()[1],
	   coo[1] = RES.Atoms()[0] - RES.Atoms()[1],
	   coo[2] = RES.Atoms()[2] - RES.Atoms()[1],
	   coo[3] = CrossProduct( x, y);

	   for( std::vector< Atom>::iterator res = RES.Atoms().begin() + 3; res != RES.Atoms().end(); ++res)
	   {
		   v.push_back( res->Position().CooTrafo( coo));
	   }

	   SafeAdd( m_Rotamers, RES.Name(), v);
   }



};


#endif /* SIDECHAIN_H_ */
