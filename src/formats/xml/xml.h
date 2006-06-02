/**********************************************************************
xml.h Declaration of XMLConversion, 
declaration and definition of XMLBaseFormat and XMLMoleculeFormat 
Copyright (C) 2005-2006 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "obconversion.h"
#include "obmolecformat.h"

#include <libxml/xmlreader.h>
#include <libxml/xmlwriter.h>

namespace OpenBabel
{

  // This macro is used in DLL builds. If it has not
  // been set in babelconfig.h, define it as nothing.
#ifndef OBCOMMON
#define OBCOMMON
#endif

  //forward declaration
  class XMLBaseFormat;

  //******************************************************
  //XMLConversion class

  /// An extended OBConversion class which includes a libxml2 reader for use with XML formats
  /// Copies an OBConversion and then extends it with a XML parser. 
  /// Instances made on the heap are deleted when the original OBConversion object is.
  class OBCOMMON XMLConversion : public OBConversion
    {
    public:
      ///Existing OBConversion instance copied
      XMLConversion(OBConversion* pConv);
		
      ///Frees reader and writer if necessary
      ~XMLConversion();

      bool SetupReader();///< opens libxml2 reader
      bool SetupWriter();///< opens libxml2 writer

      ///Parses the input xml stream and sends each element to the format's callback routines
      bool ReadXML(XMLBaseFormat* pFormat, OBBase* pOb);

      ///Read and discard XML text up to the next occurrence of the tag e.g."/molecule>"
      ///This is left as the current node. Returns 1 on success, 0 if not found, -1 if failed.
      int SkipXML(const char* ctag);

      typedef std::map<std::string, XMLBaseFormat*> NsMapType;

      ///This static function returns a reference to the map
      ///Avoids "static initialization order fiasco"
      static NsMapType& Namespaces()
        {
          static NsMapType* nsm = NULL;
          if (!nsm)
            nsm = new NsMapType;
          return *nsm;
        };

      static void RegisterXMLFormat(XMLBaseFormat* pFormat,
                                    bool IsDefault=false, const char* uri=NULL);

      ///Returns the extended OBConversion class, making it if necessary
      static XMLConversion* GetDerived(OBConversion* pConv, bool ForReading=true);

      ///Because OBConversion::Convert is still using the unextended OBConversion object
      ///we need to obtain the conversion paramters from it when requested
      bool IsLast()
        { return _pConv->IsLast(); }
      int GetOutputIndex()
        { return  _pConv->GetOutputIndex(); }


      xmlTextReaderPtr GetReader() const
        { return _reader;	};

      xmlTextWriterPtr GetWriter() const
        { return _writer;	};

      void OutputToStream()
        {
          int ret=xmlOutputBufferFlush(_buf);
        }

      static XMLBaseFormat* GetDefaultXMLClass() //TODO make dependent on object type
        { return _pDefault;};

      void LookForNamespace()
        { _LookingForNamespace = true; };

      ///Static callback functions for xmlReaderForIO()
      static int ReadStream(void * context, char * buffer, int len);
      static int WriteStream(void * context, const char * buffer, int len);
      //static int CloseStream(void* context);

      std::string GetAttribute(const char* attrname);

      ///Sets value to element content. Returns false if there is no content. 
      std::string GetContent();

      ///Sets value to element content as an integer. Returns false if there is no content. 
      bool    GetContentInt(int& value);

      ///Sets value to element content as an double. Returns false if there is no content. 
      bool GetContentDouble(double& value);

    private:
      static XMLBaseFormat* _pDefault;
      OBConversion* _pConv;
      std::streampos  _requestedpos, _lastpos;  
      xmlTextReaderPtr _reader;
      xmlTextWriterPtr _writer;
      xmlOutputBufferPtr _buf;
      //	xmlBufferPtr _buf;
      bool _LookingForNamespace;
    public:	
      bool _SkipNextRead;
    };

  //*************************************************
  /// Abstract class containing common functionality for XML formats.
  class OBCOMMON XMLBaseFormat : public OBFormat
    {
    protected:
      XMLConversion* _pxmlConv;
	
      //formating for output
      std::string _prefix;
      int baseindent, ind;
      std::string nsdecl;
      int _embedlevel;

    public:
      virtual const char* NamespaceURI()const=0;
      virtual bool DoElement(const std::string& ElName){return false;};
      virtual bool EndElement(const std::string& ElName){return false;};
      /// The tag at the end of the chemical object e.g. "/molecule>"
      virtual const char* EndTag(){return ">";};
	
    protected:
      xmlTextReaderPtr reader() const
        {
          return _pxmlConv->GetReader();
        }

      xmlTextWriterPtr writer() const
        {
          return _pxmlConv->GetWriter();
        }
	
      void OutputToStream()
        {
          _pxmlConv->OutputToStream();
        }
	
      ///Skip past first n objects in input stream (or current one with n=0)
      /// Returns 1 on success, -1 on error and 0 if not implemented 
      virtual int SkipObjects(int n, OBConversion* pConv)
        {
          //don't implement on base class
          if(*EndTag()=='>')
            return 0;

          //Set up XMLConversion class with reader 
          _pxmlConv = XMLConversion::GetDerived(pConv,true);
          if(!_pxmlConv)
            return -1;

          //always find the end of at least 1 object
          if(n==0)++n;
		
          //Skip n objects, returning -1 if not successful
          int i;
          for(i=0; i<n; ++i)
            if(_pxmlConv->SkipXML(EndTag())!=1)
              return -1;
		
          return 1;       
        }

    };

  //*************************************************
  ///Abstract class containing common functionality for XML formats which represent molecules
  class OBCOMMON XMLMoleculeFormat : public XMLBaseFormat
    {
    protected:
      OBMol* _pmol;

    public:
      virtual bool ReadChemObject(OBConversion* pConv)
        {
          return OBMoleculeFormat::ReadChemObjectImpl(pConv, this);
        };

      virtual bool WriteChemObject(OBConversion* pConv)
        {
          return OBMoleculeFormat::WriteChemObjectImpl(pConv, this);
        };

      virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv)
        {
          _pmol = dynamic_cast<OBMol*>(pOb);
          if(!_pmol)
            return false;

          _pxmlConv = XMLConversion::GetDerived(pConv,true);
          if(!_pxmlConv)
            return false;
          _embedlevel = -1;
          return _pxmlConv->ReadXML(this,pOb);
        };

      const std::type_info& GetType()
        {
          return typeid(OBMol*);
        };

    };


}//namespace

//! \file
//! \brief Declaration of XMLConversion, 
//!  declaration and definition of XMLBaseFormat and XMLMoleculeFormat 
