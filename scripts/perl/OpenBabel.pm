# This file was created automatically by SWIG 1.3.30.
# Don't modify this file, modify the SWIG interface instead.
package Chemistry::OpenBabel;
require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
package Chemistry::OpenBabelc;
bootstrap Chemistry::OpenBabel;
package Chemistry::OpenBabel;
@EXPORT = qw( ); sub dl_load_flags { 0x01 }

# ---------- BASE METHODS -------------

package Chemistry::OpenBabel;

sub TIEHASH {
    my ($classname,$obj) = @_;
    return bless $obj, $classname;
}

sub CLEAR { }

sub FIRSTKEY { }

sub NEXTKEY { }

sub FETCH {
    my ($self,$field) = @_;
    my $member_func = "swig_${field}_get";
    $self->$member_func();
}

sub STORE {
    my ($self,$field,$newval) = @_;
    my $member_func = "swig_${field}_set";
    $self->$member_func($newval);
}

sub this {
    my $ptr = shift;
    return tied(%$ptr);
}


# ------- FUNCTION WRAPPERS --------

package Chemistry::OpenBabel;

*OpenDatafile = *Chemistry::OpenBabelc::OpenDatafile;
*DoubleMultiply = *Chemistry::OpenBabelc::DoubleMultiply;
*DoubleAdd = *Chemistry::OpenBabelc::DoubleAdd;
*DoubleModulus = *Chemistry::OpenBabelc::DoubleModulus;
*__lshift__ = *Chemistry::OpenBabelc::__lshift__;
*__add__ = *Chemistry::OpenBabelc::__add__;
*__sub__ = *Chemistry::OpenBabelc::__sub__;
*__div__ = *Chemistry::OpenBabelc::__div__;
*__mul__ = *Chemistry::OpenBabelc::__mul__;
*dot = *Chemistry::OpenBabelc::dot;
*cross = *Chemistry::OpenBabelc::cross;
*vectorAngle = *Chemistry::OpenBabelc::vectorAngle;
*CalcTorsionAngle = *Chemistry::OpenBabelc::CalcTorsionAngle;
*Point2Plane = *Chemistry::OpenBabelc::Point2Plane;
*Trim = *Chemistry::OpenBabelc::Trim;
*tokenize = *Chemistry::OpenBabelc::tokenize;
*ThrowError = *Chemistry::OpenBabelc::ThrowError;
*CartesianToInternal = *Chemistry::OpenBabelc::CartesianToInternal;
*InternalToCartesian = *Chemistry::OpenBabelc::InternalToCartesian;
*NewExtension = *Chemistry::OpenBabelc::NewExtension;
*get_rmat = *Chemistry::OpenBabelc::get_rmat;
*ob_make_rmat = *Chemistry::OpenBabelc::ob_make_rmat;
*qtrfit = *Chemistry::OpenBabelc::qtrfit;
*superimpose = *Chemistry::OpenBabelc::superimpose;
*CompareRingSize = *Chemistry::OpenBabelc::CompareRingSize;
*SmartsLexReplace = *Chemistry::OpenBabelc::SmartsLexReplace;

############# Class : Chemistry::OpenBabel::vectorInt ##############

package Chemistry::OpenBabel::vectorInt;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorInt(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorInt_size;
*empty = *Chemistry::OpenBabelc::vectorInt_empty;
*clear = *Chemistry::OpenBabelc::vectorInt_clear;
*push = *Chemistry::OpenBabelc::vectorInt_push;
*pop = *Chemistry::OpenBabelc::vectorInt_pop;
*get = *Chemistry::OpenBabelc::vectorInt_get;
*set = *Chemistry::OpenBabelc::vectorInt_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorInt($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vvInt ##############

package Chemistry::OpenBabel::vvInt;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vvInt(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vvInt_size;
*empty = *Chemistry::OpenBabelc::vvInt_empty;
*clear = *Chemistry::OpenBabelc::vvInt_clear;
*push = *Chemistry::OpenBabelc::vvInt_push;
*pop = *Chemistry::OpenBabelc::vvInt_pop;
*get = *Chemistry::OpenBabelc::vvInt_get;
*set = *Chemistry::OpenBabelc::vvInt_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vvInt($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorDouble ##############

package Chemistry::OpenBabel::vectorDouble;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorDouble(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorDouble_size;
*empty = *Chemistry::OpenBabelc::vectorDouble_empty;
*clear = *Chemistry::OpenBabelc::vectorDouble_clear;
*push = *Chemistry::OpenBabelc::vectorDouble_push;
*pop = *Chemistry::OpenBabelc::vectorDouble_pop;
*get = *Chemistry::OpenBabelc::vectorDouble_get;
*set = *Chemistry::OpenBabelc::vectorDouble_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorDouble($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vVector3 ##############

package Chemistry::OpenBabel::vVector3;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vVector3(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vVector3_size;
*empty = *Chemistry::OpenBabelc::vVector3_empty;
*clear = *Chemistry::OpenBabelc::vVector3_clear;
*push = *Chemistry::OpenBabelc::vVector3_push;
*pop = *Chemistry::OpenBabelc::vVector3_pop;
*get = *Chemistry::OpenBabelc::vVector3_get;
*set = *Chemistry::OpenBabelc::vVector3_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vVector3($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorMol ##############

package Chemistry::OpenBabel::vectorMol;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorMol(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorMol_size;
*empty = *Chemistry::OpenBabelc::vectorMol_empty;
*clear = *Chemistry::OpenBabelc::vectorMol_clear;
*push = *Chemistry::OpenBabelc::vectorMol_push;
*pop = *Chemistry::OpenBabelc::vectorMol_pop;
*get = *Chemistry::OpenBabelc::vectorMol_get;
*set = *Chemistry::OpenBabelc::vectorMol_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorMol($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorBond ##############

package Chemistry::OpenBabel::vectorBond;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorBond(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorBond_size;
*empty = *Chemistry::OpenBabelc::vectorBond_empty;
*clear = *Chemistry::OpenBabelc::vectorBond_clear;
*push = *Chemistry::OpenBabelc::vectorBond_push;
*pop = *Chemistry::OpenBabelc::vectorBond_pop;
*get = *Chemistry::OpenBabelc::vectorBond_get;
*set = *Chemistry::OpenBabelc::vectorBond_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorBond($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorResidue ##############

package Chemistry::OpenBabel::vectorResidue;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorResidue(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorResidue_size;
*empty = *Chemistry::OpenBabelc::vectorResidue_empty;
*clear = *Chemistry::OpenBabelc::vectorResidue_clear;
*push = *Chemistry::OpenBabelc::vectorResidue_push;
*pop = *Chemistry::OpenBabelc::vectorResidue_pop;
*get = *Chemistry::OpenBabelc::vectorResidue_get;
*set = *Chemistry::OpenBabelc::vectorResidue_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorResidue($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorRing ##############

package Chemistry::OpenBabel::vectorRing;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorRing(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorRing_size;
*empty = *Chemistry::OpenBabelc::vectorRing_empty;
*clear = *Chemistry::OpenBabelc::vectorRing_clear;
*push = *Chemistry::OpenBabelc::vectorRing_push;
*pop = *Chemistry::OpenBabelc::vectorRing_pop;
*get = *Chemistry::OpenBabelc::vectorRing_get;
*set = *Chemistry::OpenBabelc::vectorRing_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorRing($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorData ##############

package Chemistry::OpenBabel::vectorData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorData(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorData_size;
*empty = *Chemistry::OpenBabelc::vectorData_empty;
*clear = *Chemistry::OpenBabelc::vectorData_clear;
*push = *Chemistry::OpenBabelc::vectorData_push;
*pop = *Chemistry::OpenBabelc::vectorData_pop;
*get = *Chemistry::OpenBabelc::vectorData_get;
*set = *Chemistry::OpenBabelc::vectorData_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBGlobalDataBase ##############

package Chemistry::OpenBabel::OBGlobalDataBase;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBGlobalDataBase(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBGlobalDataBase($self);
        delete $OWNER{$self};
    }
}

*Init = *Chemistry::OpenBabelc::OBGlobalDataBase_Init;
*GetSize = *Chemistry::OpenBabelc::OBGlobalDataBase_GetSize;
*SetReadDirectory = *Chemistry::OpenBabelc::OBGlobalDataBase_SetReadDirectory;
*SetEnvironmentVariable = *Chemistry::OpenBabelc::OBGlobalDataBase_SetEnvironmentVariable;
*ParseLine = *Chemistry::OpenBabelc::OBGlobalDataBase_ParseLine;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBElement ##############

package Chemistry::OpenBabel::OBElement;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBElement(@_);
    bless $self, $pkg if defined($self);
}

*GetAtomicNum = *Chemistry::OpenBabelc::OBElement_GetAtomicNum;
*GetSymbol = *Chemistry::OpenBabelc::OBElement_GetSymbol;
*GetCovalentRad = *Chemistry::OpenBabelc::OBElement_GetCovalentRad;
*GetVdwRad = *Chemistry::OpenBabelc::OBElement_GetVdwRad;
*GetMass = *Chemistry::OpenBabelc::OBElement_GetMass;
*GetMaxBonds = *Chemistry::OpenBabelc::OBElement_GetMaxBonds;
*GetElectroNeg = *Chemistry::OpenBabelc::OBElement_GetElectroNeg;
*GetIonization = *Chemistry::OpenBabelc::OBElement_GetIonization;
*GetElectronAffinity = *Chemistry::OpenBabelc::OBElement_GetElectronAffinity;
*GetName = *Chemistry::OpenBabelc::OBElement_GetName;
*GetRed = *Chemistry::OpenBabelc::OBElement_GetRed;
*GetGreen = *Chemistry::OpenBabelc::OBElement_GetGreen;
*GetBlue = *Chemistry::OpenBabelc::OBElement_GetBlue;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBElement($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBElementTable ##############

package Chemistry::OpenBabel::OBElementTable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGlobalDataBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBElementTable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBElementTable($self);
        delete $OWNER{$self};
    }
}

*ParseLine = *Chemistry::OpenBabelc::OBElementTable_ParseLine;
*GetNumberOfElements = *Chemistry::OpenBabelc::OBElementTable_GetNumberOfElements;
*GetSize = *Chemistry::OpenBabelc::OBElementTable_GetSize;
*GetAtomicNum = *Chemistry::OpenBabelc::OBElementTable_GetAtomicNum;
*GetSymbol = *Chemistry::OpenBabelc::OBElementTable_GetSymbol;
*GetVdwRad = *Chemistry::OpenBabelc::OBElementTable_GetVdwRad;
*GetCovalentRad = *Chemistry::OpenBabelc::OBElementTable_GetCovalentRad;
*GetMass = *Chemistry::OpenBabelc::OBElementTable_GetMass;
*CorrectedBondRad = *Chemistry::OpenBabelc::OBElementTable_CorrectedBondRad;
*CorrectedVdwRad = *Chemistry::OpenBabelc::OBElementTable_CorrectedVdwRad;
*GetMaxBonds = *Chemistry::OpenBabelc::OBElementTable_GetMaxBonds;
*GetElectroNeg = *Chemistry::OpenBabelc::OBElementTable_GetElectroNeg;
*GetIonization = *Chemistry::OpenBabelc::OBElementTable_GetIonization;
*GetElectronAffinity = *Chemistry::OpenBabelc::OBElementTable_GetElectronAffinity;
*GetRGB = *Chemistry::OpenBabelc::OBElementTable_GetRGB;
*GetName = *Chemistry::OpenBabelc::OBElementTable_GetName;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBIsotopeTable ##############

package Chemistry::OpenBabel::OBIsotopeTable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGlobalDataBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBIsotopeTable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBIsotopeTable($self);
        delete $OWNER{$self};
    }
}

*GetSize = *Chemistry::OpenBabelc::OBIsotopeTable_GetSize;
*ParseLine = *Chemistry::OpenBabelc::OBIsotopeTable_ParseLine;
*GetExactMass = *Chemistry::OpenBabelc::OBIsotopeTable_GetExactMass;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBTypeTable ##############

package Chemistry::OpenBabel::OBTypeTable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGlobalDataBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBTypeTable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBTypeTable($self);
        delete $OWNER{$self};
    }
}

*ParseLine = *Chemistry::OpenBabelc::OBTypeTable_ParseLine;
*GetSize = *Chemistry::OpenBabelc::OBTypeTable_GetSize;
*SetFromType = *Chemistry::OpenBabelc::OBTypeTable_SetFromType;
*SetToType = *Chemistry::OpenBabelc::OBTypeTable_SetToType;
*Translate = *Chemistry::OpenBabelc::OBTypeTable_Translate;
*GetFromType = *Chemistry::OpenBabelc::OBTypeTable_GetFromType;
*GetToType = *Chemistry::OpenBabelc::OBTypeTable_GetToType;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBResidueData ##############

package Chemistry::OpenBabel::OBResidueData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGlobalDataBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBResidueData(@_);
    bless $self, $pkg if defined($self);
}

*ParseLine = *Chemistry::OpenBabelc::OBResidueData_ParseLine;
*GetSize = *Chemistry::OpenBabelc::OBResidueData_GetSize;
*SetResName = *Chemistry::OpenBabelc::OBResidueData_SetResName;
*LookupBO = *Chemistry::OpenBabelc::OBResidueData_LookupBO;
*LookupType = *Chemistry::OpenBabelc::OBResidueData_LookupType;
*AssignBonds = *Chemistry::OpenBabelc::OBResidueData_AssignBonds;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBResidueData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBStopwatch ##############

package Chemistry::OpenBabel::OBStopwatch;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Start = *Chemistry::OpenBabelc::OBStopwatch_Start;
*Lap = *Chemistry::OpenBabelc::OBStopwatch_Lap;
*Elapsed = *Chemistry::OpenBabelc::OBStopwatch_Elapsed;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBStopwatch(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBStopwatch($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSqrtTbl ##############

package Chemistry::OpenBabel::OBSqrtTbl;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSqrtTbl(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSqrtTbl($self);
        delete $OWNER{$self};
    }
}

*Sqrt = *Chemistry::OpenBabelc::OBSqrtTbl_Sqrt;
*Init = *Chemistry::OpenBabelc::OBSqrtTbl_Init;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::DoubleType ##############

package Chemistry::OpenBabel::DoubleType;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig_hi_get = *Chemistry::OpenBabelc::DoubleType_hi_get;
*swig_hi_set = *Chemistry::OpenBabelc::DoubleType_hi_set;
*swig_lo_get = *Chemistry::OpenBabelc::DoubleType_lo_get;
*swig_lo_set = *Chemistry::OpenBabelc::DoubleType_lo_set;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_DoubleType(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_DoubleType($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRandom ##############

package Chemistry::OpenBabel::OBRandom;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRandom(@_);
    bless $self, $pkg if defined($self);
}

*Seed = *Chemistry::OpenBabelc::OBRandom_Seed;
*TimeSeed = *Chemistry::OpenBabelc::OBRandom_TimeSeed;
*NextInt = *Chemistry::OpenBabelc::OBRandom_NextInt;
*NextFloat = *Chemistry::OpenBabelc::OBRandom_NextFloat;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRandom($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vector3 ##############

package Chemistry::OpenBabel::vector3;
use overload
    "!=" => sub { $_[0]->__ne__($_[1])},
    "==" => sub { $_[0]->__eq__($_[1])},
    "fallback" => 1;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vector3(@_);
    bless $self, $pkg if defined($self);
}

*Set = *Chemistry::OpenBabelc::vector3_Set;
*SetX = *Chemistry::OpenBabelc::vector3_SetX;
*SetY = *Chemistry::OpenBabelc::vector3_SetY;
*SetZ = *Chemistry::OpenBabelc::vector3_SetZ;
*Get = *Chemistry::OpenBabelc::vector3_Get;
*AsArray = *Chemistry::OpenBabelc::vector3_AsArray;
*randomUnitVector = *Chemistry::OpenBabelc::vector3_randomUnitVector;
*normalize = *Chemistry::OpenBabelc::vector3_normalize;
*CanBeNormalized = *Chemistry::OpenBabelc::vector3_CanBeNormalized;
*length_2 = *Chemistry::OpenBabelc::vector3_length_2;
*length = *Chemistry::OpenBabelc::vector3_length;
*x = *Chemistry::OpenBabelc::vector3_x;
*y = *Chemistry::OpenBabelc::vector3_y;
*z = *Chemistry::OpenBabelc::vector3_z;
*__eq__ = *Chemistry::OpenBabelc::vector3___eq__;
*__ne__ = *Chemistry::OpenBabelc::vector3___ne__;
*IsApprox = *Chemistry::OpenBabelc::vector3_IsApprox;
*distSq = *Chemistry::OpenBabelc::vector3_distSq;
*createOrthoVector = *Chemistry::OpenBabelc::vector3_createOrthoVector;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vector3($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBGenericData ##############

package Chemistry::OpenBabel::OBGenericData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBGenericData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBGenericData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBGenericData($self);
        delete $OWNER{$self};
    }
}

*SetAttribute = *Chemistry::OpenBabelc::OBGenericData_SetAttribute;
*GetAttribute = *Chemistry::OpenBabelc::OBGenericData_GetAttribute;
*GetDataType = *Chemistry::OpenBabelc::OBGenericData_GetDataType;
*GetValue = *Chemistry::OpenBabelc::OBGenericData_GetValue;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBCommentData ##############

package Chemistry::OpenBabel::OBCommentData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBCommentData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBCommentData_Clone;
*SetData = *Chemistry::OpenBabelc::OBCommentData_SetData;
*GetData = *Chemistry::OpenBabelc::OBCommentData_GetData;
*GetValue = *Chemistry::OpenBabelc::OBCommentData_GetValue;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBCommentData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBExternalBond ##############

package Chemistry::OpenBabel::OBExternalBond;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBExternalBond(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBExternalBond($self);
        delete $OWNER{$self};
    }
}

*GetIdx = *Chemistry::OpenBabelc::OBExternalBond_GetIdx;
*GetAtom = *Chemistry::OpenBabelc::OBExternalBond_GetAtom;
*GetBond = *Chemistry::OpenBabelc::OBExternalBond_GetBond;
*SetIdx = *Chemistry::OpenBabelc::OBExternalBond_SetIdx;
*SetAtom = *Chemistry::OpenBabelc::OBExternalBond_SetAtom;
*SetBond = *Chemistry::OpenBabelc::OBExternalBond_SetBond;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBExternalBondData ##############

package Chemistry::OpenBabel::OBExternalBondData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBExternalBondData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBExternalBondData_Clone;
*SetData = *Chemistry::OpenBabelc::OBExternalBondData_SetData;
*GetData = *Chemistry::OpenBabelc::OBExternalBondData_GetData;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBExternalBondData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBPairData ##############

package Chemistry::OpenBabel::OBPairData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBPairData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBPairData_Clone;
*SetValue = *Chemistry::OpenBabelc::OBPairData_SetValue;
*GetValue = *Chemistry::OpenBabelc::OBPairData_GetValue;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBPairData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSetData ##############

package Chemistry::OpenBabel::OBSetData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSetData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBSetData_Clone;
*AddData = *Chemistry::OpenBabelc::OBSetData_AddData;
*SetData = *Chemistry::OpenBabelc::OBSetData_SetData;
*GetData = *Chemistry::OpenBabelc::OBSetData_GetData;
*GetBegin = *Chemistry::OpenBabelc::OBSetData_GetBegin;
*GetEnd = *Chemistry::OpenBabelc::OBSetData_GetEnd;
*DeleteData = *Chemistry::OpenBabelc::OBSetData_DeleteData;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSetData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBVirtualBond ##############

package Chemistry::OpenBabel::OBVirtualBond;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Clone = *Chemistry::OpenBabelc::OBVirtualBond_Clone;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBVirtualBond(@_);
    bless $self, $pkg if defined($self);
}

*GetBgn = *Chemistry::OpenBabelc::OBVirtualBond_GetBgn;
*GetEnd = *Chemistry::OpenBabelc::OBVirtualBond_GetEnd;
*GetOrder = *Chemistry::OpenBabelc::OBVirtualBond_GetOrder;
*GetStereo = *Chemistry::OpenBabelc::OBVirtualBond_GetStereo;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBVirtualBond($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRingData ##############

package Chemistry::OpenBabel::OBRingData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRingData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBRingData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRingData($self);
        delete $OWNER{$self};
    }
}

*SetData = *Chemistry::OpenBabelc::OBRingData_SetData;
*PushBack = *Chemistry::OpenBabelc::OBRingData_PushBack;
*GetData = *Chemistry::OpenBabelc::OBRingData_GetData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBUnitCell ##############

package Chemistry::OpenBabel::OBUnitCell;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Undefined = *Chemistry::OpenBabelc::OBUnitCell_Undefined;
*Triclinic = *Chemistry::OpenBabelc::OBUnitCell_Triclinic;
*Monoclinic = *Chemistry::OpenBabelc::OBUnitCell_Monoclinic;
*Orthorhombic = *Chemistry::OpenBabelc::OBUnitCell_Orthorhombic;
*Tetragonal = *Chemistry::OpenBabelc::OBUnitCell_Tetragonal;
*Rhombohedral = *Chemistry::OpenBabelc::OBUnitCell_Rhombohedral;
*Hexagonal = *Chemistry::OpenBabelc::OBUnitCell_Hexagonal;
*Cubic = *Chemistry::OpenBabelc::OBUnitCell_Cubic;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBUnitCell(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBUnitCell_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBUnitCell($self);
        delete $OWNER{$self};
    }
}

*SetData = *Chemistry::OpenBabelc::OBUnitCell_SetData;
*SetOffset = *Chemistry::OpenBabelc::OBUnitCell_SetOffset;
*SetSpaceGroup = *Chemistry::OpenBabelc::OBUnitCell_SetSpaceGroup;
*SetLatticeType = *Chemistry::OpenBabelc::OBUnitCell_SetLatticeType;
*GetA = *Chemistry::OpenBabelc::OBUnitCell_GetA;
*GetB = *Chemistry::OpenBabelc::OBUnitCell_GetB;
*GetC = *Chemistry::OpenBabelc::OBUnitCell_GetC;
*GetAlpha = *Chemistry::OpenBabelc::OBUnitCell_GetAlpha;
*GetBeta = *Chemistry::OpenBabelc::OBUnitCell_GetBeta;
*GetGamma = *Chemistry::OpenBabelc::OBUnitCell_GetGamma;
*GetOffset = *Chemistry::OpenBabelc::OBUnitCell_GetOffset;
*GetSpaceGroup = *Chemistry::OpenBabelc::OBUnitCell_GetSpaceGroup;
*GetLatticeType = *Chemistry::OpenBabelc::OBUnitCell_GetLatticeType;
*GetCellVectors = *Chemistry::OpenBabelc::OBUnitCell_GetCellVectors;
*GetCellMatrix = *Chemistry::OpenBabelc::OBUnitCell_GetCellMatrix;
*GetOrthoMatrix = *Chemistry::OpenBabelc::OBUnitCell_GetOrthoMatrix;
*GetFractionalMatrix = *Chemistry::OpenBabelc::OBUnitCell_GetFractionalMatrix;
*GetSpaceGroupNumber = *Chemistry::OpenBabelc::OBUnitCell_GetSpaceGroupNumber;
*GetCellVolume = *Chemistry::OpenBabelc::OBUnitCell_GetCellVolume;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBConformerData ##############

package Chemistry::OpenBabel::OBConformerData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBConformerData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBConformerData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBConformerData($self);
        delete $OWNER{$self};
    }
}

*SetDimension = *Chemistry::OpenBabelc::OBConformerData_SetDimension;
*SetEnergies = *Chemistry::OpenBabelc::OBConformerData_SetEnergies;
*SetForces = *Chemistry::OpenBabelc::OBConformerData_SetForces;
*SetVelocities = *Chemistry::OpenBabelc::OBConformerData_SetVelocities;
*SetDisplacements = *Chemistry::OpenBabelc::OBConformerData_SetDisplacements;
*SetData = *Chemistry::OpenBabelc::OBConformerData_SetData;
*GetDimension = *Chemistry::OpenBabelc::OBConformerData_GetDimension;
*GetEnergies = *Chemistry::OpenBabelc::OBConformerData_GetEnergies;
*GetForces = *Chemistry::OpenBabelc::OBConformerData_GetForces;
*GetVelocities = *Chemistry::OpenBabelc::OBConformerData_GetVelocities;
*GetDisplacements = *Chemistry::OpenBabelc::OBConformerData_GetDisplacements;
*GetData = *Chemistry::OpenBabelc::OBConformerData_GetData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSymmetryData ##############

package Chemistry::OpenBabel::OBSymmetryData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSymmetryData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBSymmetryData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSymmetryData($self);
        delete $OWNER{$self};
    }
}

*SetData = *Chemistry::OpenBabelc::OBSymmetryData_SetData;
*SetPointGroup = *Chemistry::OpenBabelc::OBSymmetryData_SetPointGroup;
*SetSpaceGroup = *Chemistry::OpenBabelc::OBSymmetryData_SetSpaceGroup;
*GetPointGroup = *Chemistry::OpenBabelc::OBSymmetryData_GetPointGroup;
*GetSpaceGroup = *Chemistry::OpenBabelc::OBSymmetryData_GetSpaceGroup;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBTorsion ##############

package Chemistry::OpenBabel::OBTorsion;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBTorsion(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBTorsion($self);
        delete $OWNER{$self};
    }
}

*Clear = *Chemistry::OpenBabelc::OBTorsion_Clear;
*Empty = *Chemistry::OpenBabelc::OBTorsion_Empty;
*AddTorsion = *Chemistry::OpenBabelc::OBTorsion_AddTorsion;
*SetAngle = *Chemistry::OpenBabelc::OBTorsion_SetAngle;
*SetData = *Chemistry::OpenBabelc::OBTorsion_SetData;
*GetAngle = *Chemistry::OpenBabelc::OBTorsion_GetAngle;
*GetBondIdx = *Chemistry::OpenBabelc::OBTorsion_GetBondIdx;
*GetSize = *Chemistry::OpenBabelc::OBTorsion_GetSize;
*GetBC = *Chemistry::OpenBabelc::OBTorsion_GetBC;
*GetADs = *Chemistry::OpenBabelc::OBTorsion_GetADs;
*IsProtonRotor = *Chemistry::OpenBabelc::OBTorsion_IsProtonRotor;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBTorsionData ##############

package Chemistry::OpenBabel::OBTorsionData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Clone = *Chemistry::OpenBabelc::OBTorsionData_Clone;
*Clear = *Chemistry::OpenBabelc::OBTorsionData_Clear;
*GetData = *Chemistry::OpenBabelc::OBTorsionData_GetData;
*GetSize = *Chemistry::OpenBabelc::OBTorsionData_GetSize;
*SetData = *Chemistry::OpenBabelc::OBTorsionData_SetData;
*FillTorsionArray = *Chemistry::OpenBabelc::OBTorsionData_FillTorsionArray;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBTorsionData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBAngle ##############

package Chemistry::OpenBabel::OBAngle;
use overload
    "==" => sub { $_[0]->__eq__($_[1])},
    "fallback" => 1;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBAngle(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBAngle($self);
        delete $OWNER{$self};
    }
}

*__eq__ = *Chemistry::OpenBabelc::OBAngle___eq__;
*Clear = *Chemistry::OpenBabelc::OBAngle_Clear;
*GetAngle = *Chemistry::OpenBabelc::OBAngle_GetAngle;
*SetAngle = *Chemistry::OpenBabelc::OBAngle_SetAngle;
*SetAtoms = *Chemistry::OpenBabelc::OBAngle_SetAtoms;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBAngleData ##############

package Chemistry::OpenBabel::OBAngleData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Clone = *Chemistry::OpenBabelc::OBAngleData_Clone;
*Clear = *Chemistry::OpenBabelc::OBAngleData_Clear;
*FillAngleArray = *Chemistry::OpenBabelc::OBAngleData_FillAngleArray;
*SetData = *Chemistry::OpenBabelc::OBAngleData_SetData;
*GetSize = *Chemistry::OpenBabelc::OBAngleData_GetSize;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBAngleData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBChiralData ##############

package Chemistry::OpenBabel::OBChiralData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*GetAtom4Refs = *Chemistry::OpenBabelc::OBChiralData_GetAtom4Refs;
*GetAtomRef = *Chemistry::OpenBabelc::OBChiralData_GetAtomRef;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBChiralData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBChiralData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBChiralData($self);
        delete $OWNER{$self};
    }
}

*Clear = *Chemistry::OpenBabelc::OBChiralData_Clear;
*SetAtom4Refs = *Chemistry::OpenBabelc::OBChiralData_SetAtom4Refs;
*AddAtomRef = *Chemistry::OpenBabelc::OBChiralData_AddAtomRef;
*GetSize = *Chemistry::OpenBabelc::OBChiralData_GetSize;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSerialNums ##############

package Chemistry::OpenBabel::OBSerialNums;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSerialNums(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBSerialNums_Clone;
*GetData = *Chemistry::OpenBabelc::OBSerialNums_GetData;
*SetData = *Chemistry::OpenBabelc::OBSerialNums_SetData;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSerialNums($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBBase ##############

package Chemistry::OpenBabel::OBBase;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBBase($self);
        delete $OWNER{$self};
    }
}

*DoTransformations = *Chemistry::OpenBabelc::OBBase_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBBase_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBBase_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBBase_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBBase_SetData;
*DataSize = *Chemistry::OpenBabelc::OBBase_DataSize;
*GetData = *Chemistry::OpenBabelc::OBBase_GetData;
*BeginData = *Chemistry::OpenBabelc::OBBase_BeginData;
*EndData = *Chemistry::OpenBabelc::OBBase_EndData;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBBase(@_);
    bless $self, $pkg if defined($self);
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBError ##############

package Chemistry::OpenBabel::OBError;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBError(@_);
    bless $self, $pkg if defined($self);
}

*message = *Chemistry::OpenBabelc::OBError_message;
*GetMethod = *Chemistry::OpenBabelc::OBError_GetMethod;
*GetError = *Chemistry::OpenBabelc::OBError_GetError;
*GetExplanation = *Chemistry::OpenBabelc::OBError_GetExplanation;
*GetPossibleCause = *Chemistry::OpenBabelc::OBError_GetPossibleCause;
*GetSuggestedRemedy = *Chemistry::OpenBabelc::OBError_GetSuggestedRemedy;
*GetLevel = *Chemistry::OpenBabelc::OBError_GetLevel;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBError($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBMessageHandler ##############

package Chemistry::OpenBabel::OBMessageHandler;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBMessageHandler(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBMessageHandler($self);
        delete $OWNER{$self};
    }
}

*ThrowError = *Chemistry::OpenBabelc::OBMessageHandler_ThrowError;
*GetMessagesOfLevel = *Chemistry::OpenBabelc::OBMessageHandler_GetMessagesOfLevel;
*StartLogging = *Chemistry::OpenBabelc::OBMessageHandler_StartLogging;
*StopLogging = *Chemistry::OpenBabelc::OBMessageHandler_StopLogging;
*SetMaxLogEntries = *Chemistry::OpenBabelc::OBMessageHandler_SetMaxLogEntries;
*GetMaxLogEntries = *Chemistry::OpenBabelc::OBMessageHandler_GetMaxLogEntries;
*ClearLog = *Chemistry::OpenBabelc::OBMessageHandler_ClearLog;
*SetOutputLevel = *Chemistry::OpenBabelc::OBMessageHandler_SetOutputLevel;
*GetOutputLevel = *Chemistry::OpenBabelc::OBMessageHandler_GetOutputLevel;
*SetOutputStream = *Chemistry::OpenBabelc::OBMessageHandler_SetOutputStream;
*GetOutputStream = *Chemistry::OpenBabelc::OBMessageHandler_GetOutputStream;
*StartErrorWrap = *Chemistry::OpenBabelc::OBMessageHandler_StartErrorWrap;
*StopErrorWrap = *Chemistry::OpenBabelc::OBMessageHandler_StopErrorWrap;
*GetErrorMessageCount = *Chemistry::OpenBabelc::OBMessageHandler_GetErrorMessageCount;
*GetWarningMessageCount = *Chemistry::OpenBabelc::OBMessageHandler_GetWarningMessageCount;
*GetInfoMessageCount = *Chemistry::OpenBabelc::OBMessageHandler_GetInfoMessageCount;
*GetAuditMessageCount = *Chemistry::OpenBabelc::OBMessageHandler_GetAuditMessageCount;
*GetDebugMessageCount = *Chemistry::OpenBabelc::OBMessageHandler_GetDebugMessageCount;
*GetMessageSummary = *Chemistry::OpenBabelc::OBMessageHandler_GetMessageSummary;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::obLogBuf ##############

package Chemistry::OpenBabel::obLogBuf;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_obLogBuf($self);
        delete $OWNER{$self};
    }
}

sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_obLogBuf(@_);
    bless $self, $pkg if defined($self);
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBFormat ##############

package Chemistry::OpenBabel::OBFormat;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*ReadMolecule = *Chemistry::OpenBabelc::OBFormat_ReadMolecule;
*ReadChemObject = *Chemistry::OpenBabelc::OBFormat_ReadChemObject;
*WriteMolecule = *Chemistry::OpenBabelc::OBFormat_WriteMolecule;
*WriteChemObject = *Chemistry::OpenBabelc::OBFormat_WriteChemObject;
*Description = *Chemistry::OpenBabelc::OBFormat_Description;
*TargetClassDescription = *Chemistry::OpenBabelc::OBFormat_TargetClassDescription;
*GetType = *Chemistry::OpenBabelc::OBFormat_GetType;
*SpecificationURL = *Chemistry::OpenBabelc::OBFormat_SpecificationURL;
*GetMIMEType = *Chemistry::OpenBabelc::OBFormat_GetMIMEType;
*Flags = *Chemistry::OpenBabelc::OBFormat_Flags;
*SkipObjects = *Chemistry::OpenBabelc::OBFormat_SkipObjects;
*MakeNewInstance = *Chemistry::OpenBabelc::OBFormat_MakeNewInstance;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBFormat($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::CharPtrLess ##############

package Chemistry::OpenBabel::CharPtrLess;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*__call__ = *Chemistry::OpenBabelc::CharPtrLess___call__;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_CharPtrLess(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_CharPtrLess($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBConversion ##############

package Chemistry::OpenBabel::OBConversion;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBConversion(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBConversion($self);
        delete $OWNER{$self};
    }
}

*RegisterFormat = *Chemistry::OpenBabelc::OBConversion_RegisterFormat;
*FindFormat = *Chemistry::OpenBabelc::OBConversion_FindFormat;
*FormatFromExt = *Chemistry::OpenBabelc::OBConversion_FormatFromExt;
*FormatFromMIME = *Chemistry::OpenBabelc::OBConversion_FormatFromMIME;
*GetNextFormat = *Chemistry::OpenBabelc::OBConversion_GetNextFormat;
*Description = *Chemistry::OpenBabelc::OBConversion_Description;
*GetInStream = *Chemistry::OpenBabelc::OBConversion_GetInStream;
*GetOutStream = *Chemistry::OpenBabelc::OBConversion_GetOutStream;
*SetInStream = *Chemistry::OpenBabelc::OBConversion_SetInStream;
*SetOutStream = *Chemistry::OpenBabelc::OBConversion_SetOutStream;
*SetInAndOutFormats = *Chemistry::OpenBabelc::OBConversion_SetInAndOutFormats;
*SetInFormat = *Chemistry::OpenBabelc::OBConversion_SetInFormat;
*SetOutFormat = *Chemistry::OpenBabelc::OBConversion_SetOutFormat;
*GetInFormat = *Chemistry::OpenBabelc::OBConversion_GetInFormat;
*GetOutFormat = *Chemistry::OpenBabelc::OBConversion_GetOutFormat;
*GetInFilename = *Chemistry::OpenBabelc::OBConversion_GetInFilename;
*GetInPos = *Chemistry::OpenBabelc::OBConversion_GetInPos;
*GetInLen = *Chemistry::OpenBabelc::OBConversion_GetInLen;
*GetTitle = *Chemistry::OpenBabelc::OBConversion_GetTitle;
*GetAuxConv = *Chemistry::OpenBabelc::OBConversion_GetAuxConv;
*SetAuxConv = *Chemistry::OpenBabelc::OBConversion_SetAuxConv;
*INOPTIONS = *Chemistry::OpenBabelc::OBConversion_INOPTIONS;
*OUTOPTIONS = *Chemistry::OpenBabelc::OBConversion_OUTOPTIONS;
*GENOPTIONS = *Chemistry::OpenBabelc::OBConversion_GENOPTIONS;
*IsOption = *Chemistry::OpenBabelc::OBConversion_IsOption;
*GetOptions = *Chemistry::OpenBabelc::OBConversion_GetOptions;
*AddOption = *Chemistry::OpenBabelc::OBConversion_AddOption;
*RemoveOption = *Chemistry::OpenBabelc::OBConversion_RemoveOption;
*SetOptions = *Chemistry::OpenBabelc::OBConversion_SetOptions;
*RegisterOptionParam = *Chemistry::OpenBabelc::OBConversion_RegisterOptionParam;
*GetOptionParams = *Chemistry::OpenBabelc::OBConversion_GetOptionParams;
*Convert = *Chemistry::OpenBabelc::OBConversion_Convert;
*FullConvert = *Chemistry::OpenBabelc::OBConversion_FullConvert;
*AddChemObject = *Chemistry::OpenBabelc::OBConversion_AddChemObject;
*GetChemObject = *Chemistry::OpenBabelc::OBConversion_GetChemObject;
*IsLast = *Chemistry::OpenBabelc::OBConversion_IsLast;
*IsFirstInput = *Chemistry::OpenBabelc::OBConversion_IsFirstInput;
*GetOutputIndex = *Chemistry::OpenBabelc::OBConversion_GetOutputIndex;
*SetOutputIndex = *Chemistry::OpenBabelc::OBConversion_SetOutputIndex;
*SetMoreFilesToCome = *Chemistry::OpenBabelc::OBConversion_SetMoreFilesToCome;
*SetOneObjectOnly = *Chemistry::OpenBabelc::OBConversion_SetOneObjectOnly;
*GetDefaultFormat = *Chemistry::OpenBabelc::OBConversion_GetDefaultFormat;
*Write = *Chemistry::OpenBabelc::OBConversion_Write;
*WriteString = *Chemistry::OpenBabelc::OBConversion_WriteString;
*WriteFile = *Chemistry::OpenBabelc::OBConversion_WriteFile;
*CloseOutFile = *Chemistry::OpenBabelc::OBConversion_CloseOutFile;
*Read = *Chemistry::OpenBabelc::OBConversion_Read;
*ReadString = *Chemistry::OpenBabelc::OBConversion_ReadString;
*ReadFile = *Chemistry::OpenBabelc::OBConversion_ReadFile;
*BatchFileName = *Chemistry::OpenBabelc::OBConversion_BatchFileName;
*IncrementedFileName = *Chemistry::OpenBabelc::OBConversion_IncrementedFileName;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBResidue ##############

package Chemistry::OpenBabel::OBResidue;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBResidue(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBResidue($self);
        delete $OWNER{$self};
    }
}

*AddAtom = *Chemistry::OpenBabelc::OBResidue_AddAtom;
*InsertAtom = *Chemistry::OpenBabelc::OBResidue_InsertAtom;
*RemoveAtom = *Chemistry::OpenBabelc::OBResidue_RemoveAtom;
*Clear = *Chemistry::OpenBabelc::OBResidue_Clear;
*SetName = *Chemistry::OpenBabelc::OBResidue_SetName;
*SetNum = *Chemistry::OpenBabelc::OBResidue_SetNum;
*SetChain = *Chemistry::OpenBabelc::OBResidue_SetChain;
*SetChainNum = *Chemistry::OpenBabelc::OBResidue_SetChainNum;
*SetIdx = *Chemistry::OpenBabelc::OBResidue_SetIdx;
*SetAtomID = *Chemistry::OpenBabelc::OBResidue_SetAtomID;
*SetHetAtom = *Chemistry::OpenBabelc::OBResidue_SetHetAtom;
*SetSerialNum = *Chemistry::OpenBabelc::OBResidue_SetSerialNum;
*GetName = *Chemistry::OpenBabelc::OBResidue_GetName;
*GetNum = *Chemistry::OpenBabelc::OBResidue_GetNum;
*GetNumAtoms = *Chemistry::OpenBabelc::OBResidue_GetNumAtoms;
*GetChain = *Chemistry::OpenBabelc::OBResidue_GetChain;
*GetChainNum = *Chemistry::OpenBabelc::OBResidue_GetChainNum;
*GetIdx = *Chemistry::OpenBabelc::OBResidue_GetIdx;
*GetResKey = *Chemistry::OpenBabelc::OBResidue_GetResKey;
*GetAtoms = *Chemistry::OpenBabelc::OBResidue_GetAtoms;
*GetBonds = *Chemistry::OpenBabelc::OBResidue_GetBonds;
*GetAtomID = *Chemistry::OpenBabelc::OBResidue_GetAtomID;
*GetSerialNum = *Chemistry::OpenBabelc::OBResidue_GetSerialNum;
*GetAminoAcidProperty = *Chemistry::OpenBabelc::OBResidue_GetAminoAcidProperty;
*GetAtomProperty = *Chemistry::OpenBabelc::OBResidue_GetAtomProperty;
*GetResidueProperty = *Chemistry::OpenBabelc::OBResidue_GetResidueProperty;
*IsHetAtom = *Chemistry::OpenBabelc::OBResidue_IsHetAtom;
*IsResidueType = *Chemistry::OpenBabelc::OBResidue_IsResidueType;
*BeginAtoms = *Chemistry::OpenBabelc::OBResidue_BeginAtoms;
*EndAtoms = *Chemistry::OpenBabelc::OBResidue_EndAtoms;
*BeginAtom = *Chemistry::OpenBabelc::OBResidue_BeginAtom;
*NextAtom = *Chemistry::OpenBabelc::OBResidue_NextAtom;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBAtom ##############

package Chemistry::OpenBabel::OBAtom;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig_Visit_get = *Chemistry::OpenBabelc::OBAtom_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBAtom_Visit_set;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBAtom(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBAtom($self);
        delete $OWNER{$self};
    }
}

*Clear = *Chemistry::OpenBabelc::OBAtom_Clear;
*SetIdx = *Chemistry::OpenBabelc::OBAtom_SetIdx;
*SetHyb = *Chemistry::OpenBabelc::OBAtom_SetHyb;
*SetAtomicNum = *Chemistry::OpenBabelc::OBAtom_SetAtomicNum;
*SetIsotope = *Chemistry::OpenBabelc::OBAtom_SetIsotope;
*SetImplicitValence = *Chemistry::OpenBabelc::OBAtom_SetImplicitValence;
*IncrementImplicitValence = *Chemistry::OpenBabelc::OBAtom_IncrementImplicitValence;
*DecrementImplicitValence = *Chemistry::OpenBabelc::OBAtom_DecrementImplicitValence;
*SetFormalCharge = *Chemistry::OpenBabelc::OBAtom_SetFormalCharge;
*SetSpinMultiplicity = *Chemistry::OpenBabelc::OBAtom_SetSpinMultiplicity;
*SetType = *Chemistry::OpenBabelc::OBAtom_SetType;
*SetPartialCharge = *Chemistry::OpenBabelc::OBAtom_SetPartialCharge;
*SetCoordPtr = *Chemistry::OpenBabelc::OBAtom_SetCoordPtr;
*SetVector = *Chemistry::OpenBabelc::OBAtom_SetVector;
*SetResidue = *Chemistry::OpenBabelc::OBAtom_SetResidue;
*SetParent = *Chemistry::OpenBabelc::OBAtom_SetParent;
*SetAromatic = *Chemistry::OpenBabelc::OBAtom_SetAromatic;
*UnsetAromatic = *Chemistry::OpenBabelc::OBAtom_UnsetAromatic;
*SetClockwiseStereo = *Chemistry::OpenBabelc::OBAtom_SetClockwiseStereo;
*SetAntiClockwiseStereo = *Chemistry::OpenBabelc::OBAtom_SetAntiClockwiseStereo;
*SetPositiveStereo = *Chemistry::OpenBabelc::OBAtom_SetPositiveStereo;
*SetNegativeStereo = *Chemistry::OpenBabelc::OBAtom_SetNegativeStereo;
*UnsetStereo = *Chemistry::OpenBabelc::OBAtom_UnsetStereo;
*SetInRing = *Chemistry::OpenBabelc::OBAtom_SetInRing;
*SetChiral = *Chemistry::OpenBabelc::OBAtom_SetChiral;
*ClearCoordPtr = *Chemistry::OpenBabelc::OBAtom_ClearCoordPtr;
*GetFormalCharge = *Chemistry::OpenBabelc::OBAtom_GetFormalCharge;
*GetAtomicNum = *Chemistry::OpenBabelc::OBAtom_GetAtomicNum;
*GetIsotope = *Chemistry::OpenBabelc::OBAtom_GetIsotope;
*GetSpinMultiplicity = *Chemistry::OpenBabelc::OBAtom_GetSpinMultiplicity;
*GetAtomicMass = *Chemistry::OpenBabelc::OBAtom_GetAtomicMass;
*GetExactMass = *Chemistry::OpenBabelc::OBAtom_GetExactMass;
*GetIdx = *Chemistry::OpenBabelc::OBAtom_GetIdx;
*GetCoordinateIdx = *Chemistry::OpenBabelc::OBAtom_GetCoordinateIdx;
*GetCIdx = *Chemistry::OpenBabelc::OBAtom_GetCIdx;
*GetValence = *Chemistry::OpenBabelc::OBAtom_GetValence;
*GetHyb = *Chemistry::OpenBabelc::OBAtom_GetHyb;
*GetImplicitValence = *Chemistry::OpenBabelc::OBAtom_GetImplicitValence;
*GetHvyValence = *Chemistry::OpenBabelc::OBAtom_GetHvyValence;
*GetHeteroValence = *Chemistry::OpenBabelc::OBAtom_GetHeteroValence;
*GetType = *Chemistry::OpenBabelc::OBAtom_GetType;
*GetX = *Chemistry::OpenBabelc::OBAtom_GetX;
*x = *Chemistry::OpenBabelc::OBAtom_x;
*GetY = *Chemistry::OpenBabelc::OBAtom_GetY;
*y = *Chemistry::OpenBabelc::OBAtom_y;
*GetZ = *Chemistry::OpenBabelc::OBAtom_GetZ;
*z = *Chemistry::OpenBabelc::OBAtom_z;
*GetCoordinate = *Chemistry::OpenBabelc::OBAtom_GetCoordinate;
*GetVector = *Chemistry::OpenBabelc::OBAtom_GetVector;
*GetPartialCharge = *Chemistry::OpenBabelc::OBAtom_GetPartialCharge;
*GetResidue = *Chemistry::OpenBabelc::OBAtom_GetResidue;
*GetParent = *Chemistry::OpenBabelc::OBAtom_GetParent;
*GetNewBondVector = *Chemistry::OpenBabelc::OBAtom_GetNewBondVector;
*GetBond = *Chemistry::OpenBabelc::OBAtom_GetBond;
*GetNextAtom = *Chemistry::OpenBabelc::OBAtom_GetNextAtom;
*BeginBonds = *Chemistry::OpenBabelc::OBAtom_BeginBonds;
*EndBonds = *Chemistry::OpenBabelc::OBAtom_EndBonds;
*BeginBond = *Chemistry::OpenBabelc::OBAtom_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBAtom_NextBond;
*BeginNbrAtom = *Chemistry::OpenBabelc::OBAtom_BeginNbrAtom;
*NextNbrAtom = *Chemistry::OpenBabelc::OBAtom_NextNbrAtom;
*GetDistance = *Chemistry::OpenBabelc::OBAtom_GetDistance;
*GetAngle = *Chemistry::OpenBabelc::OBAtom_GetAngle;
*NewResidue = *Chemistry::OpenBabelc::OBAtom_NewResidue;
*DeleteResidue = *Chemistry::OpenBabelc::OBAtom_DeleteResidue;
*AddBond = *Chemistry::OpenBabelc::OBAtom_AddBond;
*InsertBond = *Chemistry::OpenBabelc::OBAtom_InsertBond;
*DeleteBond = *Chemistry::OpenBabelc::OBAtom_DeleteBond;
*ClearBond = *Chemistry::OpenBabelc::OBAtom_ClearBond;
*CountFreeOxygens = *Chemistry::OpenBabelc::OBAtom_CountFreeOxygens;
*ImplicitHydrogenCount = *Chemistry::OpenBabelc::OBAtom_ImplicitHydrogenCount;
*ExplicitHydrogenCount = *Chemistry::OpenBabelc::OBAtom_ExplicitHydrogenCount;
*MemberOfRingCount = *Chemistry::OpenBabelc::OBAtom_MemberOfRingCount;
*MemberOfRingSize = *Chemistry::OpenBabelc::OBAtom_MemberOfRingSize;
*CountRingBonds = *Chemistry::OpenBabelc::OBAtom_CountRingBonds;
*SmallestBondAngle = *Chemistry::OpenBabelc::OBAtom_SmallestBondAngle;
*AverageBondAngle = *Chemistry::OpenBabelc::OBAtom_AverageBondAngle;
*BOSum = *Chemistry::OpenBabelc::OBAtom_BOSum;
*KBOSum = *Chemistry::OpenBabelc::OBAtom_KBOSum;
*HtoMethyl = *Chemistry::OpenBabelc::OBAtom_HtoMethyl;
*SetHybAndGeom = *Chemistry::OpenBabelc::OBAtom_SetHybAndGeom;
*ForceNoH = *Chemistry::OpenBabelc::OBAtom_ForceNoH;
*HasNoHForced = *Chemistry::OpenBabelc::OBAtom_HasNoHForced;
*HasResidue = *Chemistry::OpenBabelc::OBAtom_HasResidue;
*IsHydrogen = *Chemistry::OpenBabelc::OBAtom_IsHydrogen;
*IsCarbon = *Chemistry::OpenBabelc::OBAtom_IsCarbon;
*IsNitrogen = *Chemistry::OpenBabelc::OBAtom_IsNitrogen;
*IsOxygen = *Chemistry::OpenBabelc::OBAtom_IsOxygen;
*IsSulfur = *Chemistry::OpenBabelc::OBAtom_IsSulfur;
*IsPhosphorus = *Chemistry::OpenBabelc::OBAtom_IsPhosphorus;
*IsAromatic = *Chemistry::OpenBabelc::OBAtom_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBAtom_IsInRing;
*IsInRingSize = *Chemistry::OpenBabelc::OBAtom_IsInRingSize;
*IsHeteroatom = *Chemistry::OpenBabelc::OBAtom_IsHeteroatom;
*IsNotCorH = *Chemistry::OpenBabelc::OBAtom_IsNotCorH;
*IsConnected = *Chemistry::OpenBabelc::OBAtom_IsConnected;
*IsOneThree = *Chemistry::OpenBabelc::OBAtom_IsOneThree;
*IsOneFour = *Chemistry::OpenBabelc::OBAtom_IsOneFour;
*IsCarboxylOxygen = *Chemistry::OpenBabelc::OBAtom_IsCarboxylOxygen;
*IsPhosphateOxygen = *Chemistry::OpenBabelc::OBAtom_IsPhosphateOxygen;
*IsSulfateOxygen = *Chemistry::OpenBabelc::OBAtom_IsSulfateOxygen;
*IsNitroOxygen = *Chemistry::OpenBabelc::OBAtom_IsNitroOxygen;
*IsAmideNitrogen = *Chemistry::OpenBabelc::OBAtom_IsAmideNitrogen;
*IsPolarHydrogen = *Chemistry::OpenBabelc::OBAtom_IsPolarHydrogen;
*IsNonPolarHydrogen = *Chemistry::OpenBabelc::OBAtom_IsNonPolarHydrogen;
*IsAromaticNOxide = *Chemistry::OpenBabelc::OBAtom_IsAromaticNOxide;
*IsChiral = *Chemistry::OpenBabelc::OBAtom_IsChiral;
*IsAxial = *Chemistry::OpenBabelc::OBAtom_IsAxial;
*IsClockwise = *Chemistry::OpenBabelc::OBAtom_IsClockwise;
*IsAntiClockwise = *Chemistry::OpenBabelc::OBAtom_IsAntiClockwise;
*IsPositiveStereo = *Chemistry::OpenBabelc::OBAtom_IsPositiveStereo;
*IsNegativeStereo = *Chemistry::OpenBabelc::OBAtom_IsNegativeStereo;
*HasChiralitySpecified = *Chemistry::OpenBabelc::OBAtom_HasChiralitySpecified;
*HasChiralVolume = *Chemistry::OpenBabelc::OBAtom_HasChiralVolume;
*IsHbondAcceptor = *Chemistry::OpenBabelc::OBAtom_IsHbondAcceptor;
*IsHbondDonor = *Chemistry::OpenBabelc::OBAtom_IsHbondDonor;
*IsHbondDonorH = *Chemistry::OpenBabelc::OBAtom_IsHbondDonorH;
*HasAlphaBetaUnsat = *Chemistry::OpenBabelc::OBAtom_HasAlphaBetaUnsat;
*HasBondOfOrder = *Chemistry::OpenBabelc::OBAtom_HasBondOfOrder;
*CountBondsOfOrder = *Chemistry::OpenBabelc::OBAtom_CountBondsOfOrder;
*HasNonSingleBond = *Chemistry::OpenBabelc::OBAtom_HasNonSingleBond;
*HasSingleBond = *Chemistry::OpenBabelc::OBAtom_HasSingleBond;
*HasDoubleBond = *Chemistry::OpenBabelc::OBAtom_HasDoubleBond;
*HasAromaticBond = *Chemistry::OpenBabelc::OBAtom_HasAromaticBond;
*MatchesSMARTS = *Chemistry::OpenBabelc::OBAtom_MatchesSMARTS;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBBond ##############

package Chemistry::OpenBabel::OBBond;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig_Visit_get = *Chemistry::OpenBabelc::OBBond_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBBond_Visit_set;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBBond(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBBond($self);
        delete $OWNER{$self};
    }
}

*SetIdx = *Chemistry::OpenBabelc::OBBond_SetIdx;
*SetBO = *Chemistry::OpenBabelc::OBBond_SetBO;
*SetBegin = *Chemistry::OpenBabelc::OBBond_SetBegin;
*SetEnd = *Chemistry::OpenBabelc::OBBond_SetEnd;
*SetParent = *Chemistry::OpenBabelc::OBBond_SetParent;
*SetLength = *Chemistry::OpenBabelc::OBBond_SetLength;
*Set = *Chemistry::OpenBabelc::OBBond_Set;
*SetKSingle = *Chemistry::OpenBabelc::OBBond_SetKSingle;
*SetKDouble = *Chemistry::OpenBabelc::OBBond_SetKDouble;
*SetKTriple = *Chemistry::OpenBabelc::OBBond_SetKTriple;
*SetAromatic = *Chemistry::OpenBabelc::OBBond_SetAromatic;
*SetHash = *Chemistry::OpenBabelc::OBBond_SetHash;
*SetWedge = *Chemistry::OpenBabelc::OBBond_SetWedge;
*SetUp = *Chemistry::OpenBabelc::OBBond_SetUp;
*SetDown = *Chemistry::OpenBabelc::OBBond_SetDown;
*SetInRing = *Chemistry::OpenBabelc::OBBond_SetInRing;
*SetClosure = *Chemistry::OpenBabelc::OBBond_SetClosure;
*UnsetHash = *Chemistry::OpenBabelc::OBBond_UnsetHash;
*UnsetWedge = *Chemistry::OpenBabelc::OBBond_UnsetWedge;
*UnsetUp = *Chemistry::OpenBabelc::OBBond_UnsetUp;
*UnsetDown = *Chemistry::OpenBabelc::OBBond_UnsetDown;
*UnsetAromatic = *Chemistry::OpenBabelc::OBBond_UnsetAromatic;
*UnsetKekule = *Chemistry::OpenBabelc::OBBond_UnsetKekule;
*GetIdx = *Chemistry::OpenBabelc::OBBond_GetIdx;
*GetBO = *Chemistry::OpenBabelc::OBBond_GetBO;
*GetBondOrder = *Chemistry::OpenBabelc::OBBond_GetBondOrder;
*GetFlags = *Chemistry::OpenBabelc::OBBond_GetFlags;
*GetBeginAtomIdx = *Chemistry::OpenBabelc::OBBond_GetBeginAtomIdx;
*GetEndAtomIdx = *Chemistry::OpenBabelc::OBBond_GetEndAtomIdx;
*GetBeginAtom = *Chemistry::OpenBabelc::OBBond_GetBeginAtom;
*GetEndAtom = *Chemistry::OpenBabelc::OBBond_GetEndAtom;
*GetNbrAtom = *Chemistry::OpenBabelc::OBBond_GetNbrAtom;
*GetParent = *Chemistry::OpenBabelc::OBBond_GetParent;
*GetEquibLength = *Chemistry::OpenBabelc::OBBond_GetEquibLength;
*GetLength = *Chemistry::OpenBabelc::OBBond_GetLength;
*GetNbrAtomIdx = *Chemistry::OpenBabelc::OBBond_GetNbrAtomIdx;
*IsAromatic = *Chemistry::OpenBabelc::OBBond_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBBond_IsInRing;
*IsRotor = *Chemistry::OpenBabelc::OBBond_IsRotor;
*IsAmide = *Chemistry::OpenBabelc::OBBond_IsAmide;
*IsPrimaryAmide = *Chemistry::OpenBabelc::OBBond_IsPrimaryAmide;
*IsSecondaryAmide = *Chemistry::OpenBabelc::OBBond_IsSecondaryAmide;
*IsEster = *Chemistry::OpenBabelc::OBBond_IsEster;
*IsCarbonyl = *Chemistry::OpenBabelc::OBBond_IsCarbonyl;
*IsSingle = *Chemistry::OpenBabelc::OBBond_IsSingle;
*IsDouble = *Chemistry::OpenBabelc::OBBond_IsDouble;
*IsTriple = *Chemistry::OpenBabelc::OBBond_IsTriple;
*IsKSingle = *Chemistry::OpenBabelc::OBBond_IsKSingle;
*IsKDouble = *Chemistry::OpenBabelc::OBBond_IsKDouble;
*IsKTriple = *Chemistry::OpenBabelc::OBBond_IsKTriple;
*IsClosure = *Chemistry::OpenBabelc::OBBond_IsClosure;
*IsUp = *Chemistry::OpenBabelc::OBBond_IsUp;
*IsDown = *Chemistry::OpenBabelc::OBBond_IsDown;
*IsWedge = *Chemistry::OpenBabelc::OBBond_IsWedge;
*IsHash = *Chemistry::OpenBabelc::OBBond_IsHash;
*IsDoubleBondGeometry = *Chemistry::OpenBabelc::OBBond_IsDoubleBondGeometry;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBMol ##############

package Chemistry::OpenBabel::OBMol;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBMol(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBMol($self);
        delete $OWNER{$self};
    }
}

*ReserveAtoms = *Chemistry::OpenBabelc::OBMol_ReserveAtoms;
*CreateAtom = *Chemistry::OpenBabelc::OBMol_CreateAtom;
*CreateBond = *Chemistry::OpenBabelc::OBMol_CreateBond;
*CreateResidue = *Chemistry::OpenBabelc::OBMol_CreateResidue;
*DestroyAtom = *Chemistry::OpenBabelc::OBMol_DestroyAtom;
*DestroyBond = *Chemistry::OpenBabelc::OBMol_DestroyBond;
*DestroyResidue = *Chemistry::OpenBabelc::OBMol_DestroyResidue;
*AddAtom = *Chemistry::OpenBabelc::OBMol_AddAtom;
*AddBond = *Chemistry::OpenBabelc::OBMol_AddBond;
*AddResidue = *Chemistry::OpenBabelc::OBMol_AddResidue;
*InsertAtom = *Chemistry::OpenBabelc::OBMol_InsertAtom;
*DeleteAtom = *Chemistry::OpenBabelc::OBMol_DeleteAtom;
*DeleteBond = *Chemistry::OpenBabelc::OBMol_DeleteBond;
*DeleteResidue = *Chemistry::OpenBabelc::OBMol_DeleteResidue;
*NewAtom = *Chemistry::OpenBabelc::OBMol_NewAtom;
*NewBond = *Chemistry::OpenBabelc::OBMol_NewBond;
*NewResidue = *Chemistry::OpenBabelc::OBMol_NewResidue;
*BeginModify = *Chemistry::OpenBabelc::OBMol_BeginModify;
*EndModify = *Chemistry::OpenBabelc::OBMol_EndModify;
*GetMod = *Chemistry::OpenBabelc::OBMol_GetMod;
*IncrementMod = *Chemistry::OpenBabelc::OBMol_IncrementMod;
*DecrementMod = *Chemistry::OpenBabelc::OBMol_DecrementMod;
*GetFlags = *Chemistry::OpenBabelc::OBMol_GetFlags;
*GetTitle = *Chemistry::OpenBabelc::OBMol_GetTitle;
*NumAtoms = *Chemistry::OpenBabelc::OBMol_NumAtoms;
*NumBonds = *Chemistry::OpenBabelc::OBMol_NumBonds;
*NumHvyAtoms = *Chemistry::OpenBabelc::OBMol_NumHvyAtoms;
*NumResidues = *Chemistry::OpenBabelc::OBMol_NumResidues;
*NumRotors = *Chemistry::OpenBabelc::OBMol_NumRotors;
*GetAtom = *Chemistry::OpenBabelc::OBMol_GetAtom;
*GetFirstAtom = *Chemistry::OpenBabelc::OBMol_GetFirstAtom;
*GetBond = *Chemistry::OpenBabelc::OBMol_GetBond;
*GetResidue = *Chemistry::OpenBabelc::OBMol_GetResidue;
*GetInternalCoord = *Chemistry::OpenBabelc::OBMol_GetInternalCoord;
*GetTorsion = *Chemistry::OpenBabelc::OBMol_GetTorsion;
*GetAngle = *Chemistry::OpenBabelc::OBMol_GetAngle;
*GetFormula = *Chemistry::OpenBabelc::OBMol_GetFormula;
*GetSpacedFormula = *Chemistry::OpenBabelc::OBMol_GetSpacedFormula;
*GetEnergy = *Chemistry::OpenBabelc::OBMol_GetEnergy;
*GetMolWt = *Chemistry::OpenBabelc::OBMol_GetMolWt;
*GetExactMass = *Chemistry::OpenBabelc::OBMol_GetExactMass;
*GetTotalCharge = *Chemistry::OpenBabelc::OBMol_GetTotalCharge;
*GetTotalSpinMultiplicity = *Chemistry::OpenBabelc::OBMol_GetTotalSpinMultiplicity;
*GetDimension = *Chemistry::OpenBabelc::OBMol_GetDimension;
*GetCoordinates = *Chemistry::OpenBabelc::OBMol_GetCoordinates;
*GetSSSR = *Chemistry::OpenBabelc::OBMol_GetSSSR;
*AutomaticFormalCharge = *Chemistry::OpenBabelc::OBMol_AutomaticFormalCharge;
*AutomaticPartialCharge = *Chemistry::OpenBabelc::OBMol_AutomaticPartialCharge;
*SetTitle = *Chemistry::OpenBabelc::OBMol_SetTitle;
*SetFormula = *Chemistry::OpenBabelc::OBMol_SetFormula;
*SetEnergy = *Chemistry::OpenBabelc::OBMol_SetEnergy;
*SetDimension = *Chemistry::OpenBabelc::OBMol_SetDimension;
*SetTotalCharge = *Chemistry::OpenBabelc::OBMol_SetTotalCharge;
*SetTotalSpinMultiplicity = *Chemistry::OpenBabelc::OBMol_SetTotalSpinMultiplicity;
*SetInternalCoord = *Chemistry::OpenBabelc::OBMol_SetInternalCoord;
*SetAutomaticFormalCharge = *Chemistry::OpenBabelc::OBMol_SetAutomaticFormalCharge;
*SetAutomaticPartialCharge = *Chemistry::OpenBabelc::OBMol_SetAutomaticPartialCharge;
*SetAromaticPerceived = *Chemistry::OpenBabelc::OBMol_SetAromaticPerceived;
*SetSSSRPerceived = *Chemistry::OpenBabelc::OBMol_SetSSSRPerceived;
*SetRingAtomsAndBondsPerceived = *Chemistry::OpenBabelc::OBMol_SetRingAtomsAndBondsPerceived;
*SetAtomTypesPerceived = *Chemistry::OpenBabelc::OBMol_SetAtomTypesPerceived;
*SetChainsPerceived = *Chemistry::OpenBabelc::OBMol_SetChainsPerceived;
*SetChiralityPerceived = *Chemistry::OpenBabelc::OBMol_SetChiralityPerceived;
*SetPartialChargesPerceived = *Chemistry::OpenBabelc::OBMol_SetPartialChargesPerceived;
*SetHybridizationPerceived = *Chemistry::OpenBabelc::OBMol_SetHybridizationPerceived;
*SetImplicitValencePerceived = *Chemistry::OpenBabelc::OBMol_SetImplicitValencePerceived;
*SetKekulePerceived = *Chemistry::OpenBabelc::OBMol_SetKekulePerceived;
*SetClosureBondsPerceived = *Chemistry::OpenBabelc::OBMol_SetClosureBondsPerceived;
*SetHydrogensAdded = *Chemistry::OpenBabelc::OBMol_SetHydrogensAdded;
*SetCorrectedForPH = *Chemistry::OpenBabelc::OBMol_SetCorrectedForPH;
*SetAromaticCorrected = *Chemistry::OpenBabelc::OBMol_SetAromaticCorrected;
*SetSpinMultiplicityAssigned = *Chemistry::OpenBabelc::OBMol_SetSpinMultiplicityAssigned;
*SetFlags = *Chemistry::OpenBabelc::OBMol_SetFlags;
*UnsetAromaticPerceived = *Chemistry::OpenBabelc::OBMol_UnsetAromaticPerceived;
*UnsetPartialChargesPerceived = *Chemistry::OpenBabelc::OBMol_UnsetPartialChargesPerceived;
*UnsetImplicitValencePerceived = *Chemistry::OpenBabelc::OBMol_UnsetImplicitValencePerceived;
*UnsetFlag = *Chemistry::OpenBabelc::OBMol_UnsetFlag;
*DoTransformations = *Chemistry::OpenBabelc::OBMol_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBMol_ClassDescription;
*Clear = *Chemistry::OpenBabelc::OBMol_Clear;
*RenumberAtoms = *Chemistry::OpenBabelc::OBMol_RenumberAtoms;
*ToInertialFrame = *Chemistry::OpenBabelc::OBMol_ToInertialFrame;
*Translate = *Chemistry::OpenBabelc::OBMol_Translate;
*Rotate = *Chemistry::OpenBabelc::OBMol_Rotate;
*Kekulize = *Chemistry::OpenBabelc::OBMol_Kekulize;
*PerceiveKekuleBonds = *Chemistry::OpenBabelc::OBMol_PerceiveKekuleBonds;
*NewPerceiveKekuleBonds = *Chemistry::OpenBabelc::OBMol_NewPerceiveKekuleBonds;
*DeleteHydrogen = *Chemistry::OpenBabelc::OBMol_DeleteHydrogen;
*DeleteHydrogens = *Chemistry::OpenBabelc::OBMol_DeleteHydrogens;
*DeleteNonPolarHydrogens = *Chemistry::OpenBabelc::OBMol_DeleteNonPolarHydrogens;
*AddHydrogens = *Chemistry::OpenBabelc::OBMol_AddHydrogens;
*AddPolarHydrogens = *Chemistry::OpenBabelc::OBMol_AddPolarHydrogens;
*StripSalts = *Chemistry::OpenBabelc::OBMol_StripSalts;
*ConvertDativeBonds = *Chemistry::OpenBabelc::OBMol_ConvertDativeBonds;
*CorrectForPH = *Chemistry::OpenBabelc::OBMol_CorrectForPH;
*AssignSpinMultiplicity = *Chemistry::OpenBabelc::OBMol_AssignSpinMultiplicity;
*Center = *Chemistry::OpenBabelc::OBMol_Center;
*SetTorsion = *Chemistry::OpenBabelc::OBMol_SetTorsion;
*FindSSSR = *Chemistry::OpenBabelc::OBMol_FindSSSR;
*FindRingAtomsAndBonds = *Chemistry::OpenBabelc::OBMol_FindRingAtomsAndBonds;
*FindChiralCenters = *Chemistry::OpenBabelc::OBMol_FindChiralCenters;
*FindChildren = *Chemistry::OpenBabelc::OBMol_FindChildren;
*FindLargestFragment = *Chemistry::OpenBabelc::OBMol_FindLargestFragment;
*ContigFragList = *Chemistry::OpenBabelc::OBMol_ContigFragList;
*Align = *Chemistry::OpenBabelc::OBMol_Align;
*ConnectTheDots = *Chemistry::OpenBabelc::OBMol_ConnectTheDots;
*PerceiveBondOrders = *Chemistry::OpenBabelc::OBMol_PerceiveBondOrders;
*FindTorsions = *Chemistry::OpenBabelc::OBMol_FindTorsions;
*GetGTDVector = *Chemistry::OpenBabelc::OBMol_GetGTDVector;
*GetGIVector = *Chemistry::OpenBabelc::OBMol_GetGIVector;
*GetGIDVector = *Chemistry::OpenBabelc::OBMol_GetGIDVector;
*Has2D = *Chemistry::OpenBabelc::OBMol_Has2D;
*Has3D = *Chemistry::OpenBabelc::OBMol_Has3D;
*HasNonZeroCoords = *Chemistry::OpenBabelc::OBMol_HasNonZeroCoords;
*HasAromaticPerceived = *Chemistry::OpenBabelc::OBMol_HasAromaticPerceived;
*HasSSSRPerceived = *Chemistry::OpenBabelc::OBMol_HasSSSRPerceived;
*HasRingAtomsAndBondsPerceived = *Chemistry::OpenBabelc::OBMol_HasRingAtomsAndBondsPerceived;
*HasAtomTypesPerceived = *Chemistry::OpenBabelc::OBMol_HasAtomTypesPerceived;
*HasChiralityPerceived = *Chemistry::OpenBabelc::OBMol_HasChiralityPerceived;
*HasPartialChargesPerceived = *Chemistry::OpenBabelc::OBMol_HasPartialChargesPerceived;
*HasHybridizationPerceived = *Chemistry::OpenBabelc::OBMol_HasHybridizationPerceived;
*HasImplicitValencePerceived = *Chemistry::OpenBabelc::OBMol_HasImplicitValencePerceived;
*HasKekulePerceived = *Chemistry::OpenBabelc::OBMol_HasKekulePerceived;
*HasClosureBondsPerceived = *Chemistry::OpenBabelc::OBMol_HasClosureBondsPerceived;
*HasChainsPerceived = *Chemistry::OpenBabelc::OBMol_HasChainsPerceived;
*HasHydrogensAdded = *Chemistry::OpenBabelc::OBMol_HasHydrogensAdded;
*HasAromaticCorrected = *Chemistry::OpenBabelc::OBMol_HasAromaticCorrected;
*IsCorrectedForPH = *Chemistry::OpenBabelc::OBMol_IsCorrectedForPH;
*HasSpinMultiplicityAssigned = *Chemistry::OpenBabelc::OBMol_HasSpinMultiplicityAssigned;
*IsChiral = *Chemistry::OpenBabelc::OBMol_IsChiral;
*Empty = *Chemistry::OpenBabelc::OBMol_Empty;
*NumConformers = *Chemistry::OpenBabelc::OBMol_NumConformers;
*SetConformers = *Chemistry::OpenBabelc::OBMol_SetConformers;
*AddConformer = *Chemistry::OpenBabelc::OBMol_AddConformer;
*SetConformer = *Chemistry::OpenBabelc::OBMol_SetConformer;
*CopyConformer = *Chemistry::OpenBabelc::OBMol_CopyConformer;
*DeleteConformer = *Chemistry::OpenBabelc::OBMol_DeleteConformer;
*GetConformer = *Chemistry::OpenBabelc::OBMol_GetConformer;
*BeginConformer = *Chemistry::OpenBabelc::OBMol_BeginConformer;
*NextConformer = *Chemistry::OpenBabelc::OBMol_NextConformer;
*GetConformers = *Chemistry::OpenBabelc::OBMol_GetConformers;
*BeginAtoms = *Chemistry::OpenBabelc::OBMol_BeginAtoms;
*EndAtoms = *Chemistry::OpenBabelc::OBMol_EndAtoms;
*BeginBonds = *Chemistry::OpenBabelc::OBMol_BeginBonds;
*EndBonds = *Chemistry::OpenBabelc::OBMol_EndBonds;
*BeginResidues = *Chemistry::OpenBabelc::OBMol_BeginResidues;
*EndResidues = *Chemistry::OpenBabelc::OBMol_EndResidues;
*BeginAtom = *Chemistry::OpenBabelc::OBMol_BeginAtom;
*NextAtom = *Chemistry::OpenBabelc::OBMol_NextAtom;
*BeginBond = *Chemistry::OpenBabelc::OBMol_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBMol_NextBond;
*BeginResidue = *Chemistry::OpenBabelc::OBMol_BeginResidue;
*NextResidue = *Chemistry::OpenBabelc::OBMol_NextResidue;
*BeginInternalCoord = *Chemistry::OpenBabelc::OBMol_BeginInternalCoord;
*NextInternalCoord = *Chemistry::OpenBabelc::OBMol_NextInternalCoord;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRTree ##############

package Chemistry::OpenBabel::OBRTree;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRTree(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRTree($self);
        delete $OWNER{$self};
    }
}

*GetAtomIdx = *Chemistry::OpenBabelc::OBRTree_GetAtomIdx;
*PathToRoot = *Chemistry::OpenBabelc::OBRTree_PathToRoot;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRing ##############

package Chemistry::OpenBabel::OBRing;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig__path_get = *Chemistry::OpenBabelc::OBRing__path_get;
*swig__path_set = *Chemistry::OpenBabelc::OBRing__path_set;
*swig__pathset_get = *Chemistry::OpenBabelc::OBRing__pathset_get;
*swig__pathset_set = *Chemistry::OpenBabelc::OBRing__pathset_set;
*findCenterAndNormal = *Chemistry::OpenBabelc::OBRing_findCenterAndNormal;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRing(@_);
    bless $self, $pkg if defined($self);
}

*Size = *Chemistry::OpenBabelc::OBRing_Size;
*PathSize = *Chemistry::OpenBabelc::OBRing_PathSize;
*IsMember = *Chemistry::OpenBabelc::OBRing_IsMember;
*IsAromatic = *Chemistry::OpenBabelc::OBRing_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBRing_IsInRing;
*SetParent = *Chemistry::OpenBabelc::OBRing_SetParent;
*GetParent = *Chemistry::OpenBabelc::OBRing_GetParent;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRing($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRingSearch ##############

package Chemistry::OpenBabel::OBRingSearch;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRingSearch(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRingSearch($self);
        delete $OWNER{$self};
    }
}

*SortRings = *Chemistry::OpenBabelc::OBRingSearch_SortRings;
*RemoveRedundant = *Chemistry::OpenBabelc::OBRingSearch_RemoveRedundant;
*AddRingFromClosure = *Chemistry::OpenBabelc::OBRingSearch_AddRingFromClosure;
*WriteRings = *Chemistry::OpenBabelc::OBRingSearch_WriteRings;
*SaveUniqueRing = *Chemistry::OpenBabelc::OBRingSearch_SaveUniqueRing;
*BeginRings = *Chemistry::OpenBabelc::OBRingSearch_BeginRings;
*EndRings = *Chemistry::OpenBabelc::OBRingSearch_EndRings;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSmartsPattern ##############

package Chemistry::OpenBabel::OBSmartsPattern;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSmartsPattern($self);
        delete $OWNER{$self};
    }
}

sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSmartsPattern(@_);
    bless $self, $pkg if defined($self);
}

*NumMatches = *Chemistry::OpenBabelc::OBSmartsPattern_NumMatches;
*NumAtoms = *Chemistry::OpenBabelc::OBSmartsPattern_NumAtoms;
*NumBonds = *Chemistry::OpenBabelc::OBSmartsPattern_NumBonds;
*GetAtomicNum = *Chemistry::OpenBabelc::OBSmartsPattern_GetAtomicNum;
*GetBond = *Chemistry::OpenBabelc::OBSmartsPattern_GetBond;
*GetCharge = *Chemistry::OpenBabelc::OBSmartsPattern_GetCharge;
*GetSMARTS = *Chemistry::OpenBabelc::OBSmartsPattern_GetSMARTS;
*GetVectorBinding = *Chemistry::OpenBabelc::OBSmartsPattern_GetVectorBinding;
*Empty = *Chemistry::OpenBabelc::OBSmartsPattern_Empty;
*IsValid = *Chemistry::OpenBabelc::OBSmartsPattern_IsValid;
*Init = *Chemistry::OpenBabelc::OBSmartsPattern_Init;
*WriteMapList = *Chemistry::OpenBabelc::OBSmartsPattern_WriteMapList;
*Match = *Chemistry::OpenBabelc::OBSmartsPattern_Match;
*RestrictedMatch = *Chemistry::OpenBabelc::OBSmartsPattern_RestrictedMatch;
*GetMapList = *Chemistry::OpenBabelc::OBSmartsPattern_GetMapList;
*GetUMapList = *Chemistry::OpenBabelc::OBSmartsPattern_GetUMapList;
*BeginMList = *Chemistry::OpenBabelc::OBSmartsPattern_BeginMList;
*EndMList = *Chemistry::OpenBabelc::OBSmartsPattern_EndMList;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSSMatch ##############

package Chemistry::OpenBabel::OBSSMatch;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSSMatch(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSSMatch($self);
        delete $OWNER{$self};
    }
}

*Match = *Chemistry::OpenBabelc::OBSSMatch_Match;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBMolAtomIter ##############

package Chemistry::OpenBabel::OBMolAtomIter;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBMolAtomIter(@_);
    bless $self, $pkg if defined($self);
}

*good = *Chemistry::OpenBabelc::OBMolAtomIter_good;
*inc = *Chemistry::OpenBabelc::OBMolAtomIter_inc;
*deref = *Chemistry::OpenBabelc::OBMolAtomIter_deref;
*__ref__ = *Chemistry::OpenBabelc::OBMolAtomIter___ref__;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBMolAtomIter($self);
        delete $OWNER{$self};
    }
}

*swig_Visit_get = *Chemistry::OpenBabelc::OBMolAtomIter_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBMolAtomIter_Visit_set;
*Clear = *Chemistry::OpenBabelc::OBMolAtomIter_Clear;
*SetIdx = *Chemistry::OpenBabelc::OBMolAtomIter_SetIdx;
*SetHyb = *Chemistry::OpenBabelc::OBMolAtomIter_SetHyb;
*SetAtomicNum = *Chemistry::OpenBabelc::OBMolAtomIter_SetAtomicNum;
*SetIsotope = *Chemistry::OpenBabelc::OBMolAtomIter_SetIsotope;
*SetImplicitValence = *Chemistry::OpenBabelc::OBMolAtomIter_SetImplicitValence;
*IncrementImplicitValence = *Chemistry::OpenBabelc::OBMolAtomIter_IncrementImplicitValence;
*DecrementImplicitValence = *Chemistry::OpenBabelc::OBMolAtomIter_DecrementImplicitValence;
*SetFormalCharge = *Chemistry::OpenBabelc::OBMolAtomIter_SetFormalCharge;
*SetSpinMultiplicity = *Chemistry::OpenBabelc::OBMolAtomIter_SetSpinMultiplicity;
*SetType = *Chemistry::OpenBabelc::OBMolAtomIter_SetType;
*SetPartialCharge = *Chemistry::OpenBabelc::OBMolAtomIter_SetPartialCharge;
*SetVector = *Chemistry::OpenBabelc::OBMolAtomIter_SetVector;
*SetCoordPtr = *Chemistry::OpenBabelc::OBMolAtomIter_SetCoordPtr;
*SetResidue = *Chemistry::OpenBabelc::OBMolAtomIter_SetResidue;
*SetParent = *Chemistry::OpenBabelc::OBMolAtomIter_SetParent;
*SetAromatic = *Chemistry::OpenBabelc::OBMolAtomIter_SetAromatic;
*UnsetAromatic = *Chemistry::OpenBabelc::OBMolAtomIter_UnsetAromatic;
*SetClockwiseStereo = *Chemistry::OpenBabelc::OBMolAtomIter_SetClockwiseStereo;
*SetAntiClockwiseStereo = *Chemistry::OpenBabelc::OBMolAtomIter_SetAntiClockwiseStereo;
*SetPositiveStereo = *Chemistry::OpenBabelc::OBMolAtomIter_SetPositiveStereo;
*SetNegativeStereo = *Chemistry::OpenBabelc::OBMolAtomIter_SetNegativeStereo;
*UnsetStereo = *Chemistry::OpenBabelc::OBMolAtomIter_UnsetStereo;
*SetInRing = *Chemistry::OpenBabelc::OBMolAtomIter_SetInRing;
*SetChiral = *Chemistry::OpenBabelc::OBMolAtomIter_SetChiral;
*ClearCoordPtr = *Chemistry::OpenBabelc::OBMolAtomIter_ClearCoordPtr;
*GetFormalCharge = *Chemistry::OpenBabelc::OBMolAtomIter_GetFormalCharge;
*GetAtomicNum = *Chemistry::OpenBabelc::OBMolAtomIter_GetAtomicNum;
*GetIsotope = *Chemistry::OpenBabelc::OBMolAtomIter_GetIsotope;
*GetSpinMultiplicity = *Chemistry::OpenBabelc::OBMolAtomIter_GetSpinMultiplicity;
*GetAtomicMass = *Chemistry::OpenBabelc::OBMolAtomIter_GetAtomicMass;
*GetExactMass = *Chemistry::OpenBabelc::OBMolAtomIter_GetExactMass;
*GetIdx = *Chemistry::OpenBabelc::OBMolAtomIter_GetIdx;
*GetCoordinateIdx = *Chemistry::OpenBabelc::OBMolAtomIter_GetCoordinateIdx;
*GetCIdx = *Chemistry::OpenBabelc::OBMolAtomIter_GetCIdx;
*GetValence = *Chemistry::OpenBabelc::OBMolAtomIter_GetValence;
*GetHyb = *Chemistry::OpenBabelc::OBMolAtomIter_GetHyb;
*GetImplicitValence = *Chemistry::OpenBabelc::OBMolAtomIter_GetImplicitValence;
*GetHvyValence = *Chemistry::OpenBabelc::OBMolAtomIter_GetHvyValence;
*GetHeteroValence = *Chemistry::OpenBabelc::OBMolAtomIter_GetHeteroValence;
*GetType = *Chemistry::OpenBabelc::OBMolAtomIter_GetType;
*GetX = *Chemistry::OpenBabelc::OBMolAtomIter_GetX;
*x = *Chemistry::OpenBabelc::OBMolAtomIter_x;
*GetY = *Chemistry::OpenBabelc::OBMolAtomIter_GetY;
*y = *Chemistry::OpenBabelc::OBMolAtomIter_y;
*GetZ = *Chemistry::OpenBabelc::OBMolAtomIter_GetZ;
*z = *Chemistry::OpenBabelc::OBMolAtomIter_z;
*GetCoordinate = *Chemistry::OpenBabelc::OBMolAtomIter_GetCoordinate;
*GetVector = *Chemistry::OpenBabelc::OBMolAtomIter_GetVector;
*GetPartialCharge = *Chemistry::OpenBabelc::OBMolAtomIter_GetPartialCharge;
*GetResidue = *Chemistry::OpenBabelc::OBMolAtomIter_GetResidue;
*GetParent = *Chemistry::OpenBabelc::OBMolAtomIter_GetParent;
*GetNewBondVector = *Chemistry::OpenBabelc::OBMolAtomIter_GetNewBondVector;
*GetBond = *Chemistry::OpenBabelc::OBMolAtomIter_GetBond;
*GetNextAtom = *Chemistry::OpenBabelc::OBMolAtomIter_GetNextAtom;
*BeginBonds = *Chemistry::OpenBabelc::OBMolAtomIter_BeginBonds;
*EndBonds = *Chemistry::OpenBabelc::OBMolAtomIter_EndBonds;
*BeginBond = *Chemistry::OpenBabelc::OBMolAtomIter_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBMolAtomIter_NextBond;
*BeginNbrAtom = *Chemistry::OpenBabelc::OBMolAtomIter_BeginNbrAtom;
*NextNbrAtom = *Chemistry::OpenBabelc::OBMolAtomIter_NextNbrAtom;
*GetDistance = *Chemistry::OpenBabelc::OBMolAtomIter_GetDistance;
*GetAngle = *Chemistry::OpenBabelc::OBMolAtomIter_GetAngle;
*NewResidue = *Chemistry::OpenBabelc::OBMolAtomIter_NewResidue;
*DeleteResidue = *Chemistry::OpenBabelc::OBMolAtomIter_DeleteResidue;
*AddBond = *Chemistry::OpenBabelc::OBMolAtomIter_AddBond;
*InsertBond = *Chemistry::OpenBabelc::OBMolAtomIter_InsertBond;
*DeleteBond = *Chemistry::OpenBabelc::OBMolAtomIter_DeleteBond;
*ClearBond = *Chemistry::OpenBabelc::OBMolAtomIter_ClearBond;
*CountFreeOxygens = *Chemistry::OpenBabelc::OBMolAtomIter_CountFreeOxygens;
*ImplicitHydrogenCount = *Chemistry::OpenBabelc::OBMolAtomIter_ImplicitHydrogenCount;
*ExplicitHydrogenCount = *Chemistry::OpenBabelc::OBMolAtomIter_ExplicitHydrogenCount;
*MemberOfRingCount = *Chemistry::OpenBabelc::OBMolAtomIter_MemberOfRingCount;
*MemberOfRingSize = *Chemistry::OpenBabelc::OBMolAtomIter_MemberOfRingSize;
*CountRingBonds = *Chemistry::OpenBabelc::OBMolAtomIter_CountRingBonds;
*SmallestBondAngle = *Chemistry::OpenBabelc::OBMolAtomIter_SmallestBondAngle;
*AverageBondAngle = *Chemistry::OpenBabelc::OBMolAtomIter_AverageBondAngle;
*BOSum = *Chemistry::OpenBabelc::OBMolAtomIter_BOSum;
*KBOSum = *Chemistry::OpenBabelc::OBMolAtomIter_KBOSum;
*HtoMethyl = *Chemistry::OpenBabelc::OBMolAtomIter_HtoMethyl;
*SetHybAndGeom = *Chemistry::OpenBabelc::OBMolAtomIter_SetHybAndGeom;
*ForceNoH = *Chemistry::OpenBabelc::OBMolAtomIter_ForceNoH;
*HasNoHForced = *Chemistry::OpenBabelc::OBMolAtomIter_HasNoHForced;
*HasResidue = *Chemistry::OpenBabelc::OBMolAtomIter_HasResidue;
*IsHydrogen = *Chemistry::OpenBabelc::OBMolAtomIter_IsHydrogen;
*IsCarbon = *Chemistry::OpenBabelc::OBMolAtomIter_IsCarbon;
*IsNitrogen = *Chemistry::OpenBabelc::OBMolAtomIter_IsNitrogen;
*IsOxygen = *Chemistry::OpenBabelc::OBMolAtomIter_IsOxygen;
*IsSulfur = *Chemistry::OpenBabelc::OBMolAtomIter_IsSulfur;
*IsPhosphorus = *Chemistry::OpenBabelc::OBMolAtomIter_IsPhosphorus;
*IsAromatic = *Chemistry::OpenBabelc::OBMolAtomIter_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBMolAtomIter_IsInRing;
*IsInRingSize = *Chemistry::OpenBabelc::OBMolAtomIter_IsInRingSize;
*IsHeteroatom = *Chemistry::OpenBabelc::OBMolAtomIter_IsHeteroatom;
*IsNotCorH = *Chemistry::OpenBabelc::OBMolAtomIter_IsNotCorH;
*IsConnected = *Chemistry::OpenBabelc::OBMolAtomIter_IsConnected;
*IsOneThree = *Chemistry::OpenBabelc::OBMolAtomIter_IsOneThree;
*IsOneFour = *Chemistry::OpenBabelc::OBMolAtomIter_IsOneFour;
*IsCarboxylOxygen = *Chemistry::OpenBabelc::OBMolAtomIter_IsCarboxylOxygen;
*IsPhosphateOxygen = *Chemistry::OpenBabelc::OBMolAtomIter_IsPhosphateOxygen;
*IsSulfateOxygen = *Chemistry::OpenBabelc::OBMolAtomIter_IsSulfateOxygen;
*IsNitroOxygen = *Chemistry::OpenBabelc::OBMolAtomIter_IsNitroOxygen;
*IsAmideNitrogen = *Chemistry::OpenBabelc::OBMolAtomIter_IsAmideNitrogen;
*IsPolarHydrogen = *Chemistry::OpenBabelc::OBMolAtomIter_IsPolarHydrogen;
*IsNonPolarHydrogen = *Chemistry::OpenBabelc::OBMolAtomIter_IsNonPolarHydrogen;
*IsAromaticNOxide = *Chemistry::OpenBabelc::OBMolAtomIter_IsAromaticNOxide;
*IsChiral = *Chemistry::OpenBabelc::OBMolAtomIter_IsChiral;
*IsAxial = *Chemistry::OpenBabelc::OBMolAtomIter_IsAxial;
*IsClockwise = *Chemistry::OpenBabelc::OBMolAtomIter_IsClockwise;
*IsAntiClockwise = *Chemistry::OpenBabelc::OBMolAtomIter_IsAntiClockwise;
*IsPositiveStereo = *Chemistry::OpenBabelc::OBMolAtomIter_IsPositiveStereo;
*IsNegativeStereo = *Chemistry::OpenBabelc::OBMolAtomIter_IsNegativeStereo;
*HasChiralitySpecified = *Chemistry::OpenBabelc::OBMolAtomIter_HasChiralitySpecified;
*HasChiralVolume = *Chemistry::OpenBabelc::OBMolAtomIter_HasChiralVolume;
*IsHbondAcceptor = *Chemistry::OpenBabelc::OBMolAtomIter_IsHbondAcceptor;
*IsHbondDonor = *Chemistry::OpenBabelc::OBMolAtomIter_IsHbondDonor;
*IsHbondDonorH = *Chemistry::OpenBabelc::OBMolAtomIter_IsHbondDonorH;
*HasAlphaBetaUnsat = *Chemistry::OpenBabelc::OBMolAtomIter_HasAlphaBetaUnsat;
*HasBondOfOrder = *Chemistry::OpenBabelc::OBMolAtomIter_HasBondOfOrder;
*CountBondsOfOrder = *Chemistry::OpenBabelc::OBMolAtomIter_CountBondsOfOrder;
*HasNonSingleBond = *Chemistry::OpenBabelc::OBMolAtomIter_HasNonSingleBond;
*HasSingleBond = *Chemistry::OpenBabelc::OBMolAtomIter_HasSingleBond;
*HasDoubleBond = *Chemistry::OpenBabelc::OBMolAtomIter_HasDoubleBond;
*HasAromaticBond = *Chemistry::OpenBabelc::OBMolAtomIter_HasAromaticBond;
*MatchesSMARTS = *Chemistry::OpenBabelc::OBMolAtomIter_MatchesSMARTS;
*DoTransformations = *Chemistry::OpenBabelc::OBMolAtomIter_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBMolAtomIter_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBMolAtomIter_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBMolAtomIter_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBMolAtomIter_SetData;
*DataSize = *Chemistry::OpenBabelc::OBMolAtomIter_DataSize;
*GetData = *Chemistry::OpenBabelc::OBMolAtomIter_GetData;
*BeginData = *Chemistry::OpenBabelc::OBMolAtomIter_BeginData;
*EndData = *Chemistry::OpenBabelc::OBMolAtomIter_EndData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBMolAtomDFSIter ##############

package Chemistry::OpenBabel::OBMolAtomDFSIter;
use overload
    "++" => sub { $_[0]->__plusplus__()},
    "fallback" => 1;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBMolAtomDFSIter(@_);
    bless $self, $pkg if defined($self);
}

*__plusplus__ = *Chemistry::OpenBabelc::OBMolAtomDFSIter___plusplus__;
*__deref__ = *Chemistry::OpenBabelc::OBMolAtomDFSIter___deref__;
*__ref__ = *Chemistry::OpenBabelc::OBMolAtomDFSIter___ref__;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBMolAtomDFSIter($self);
        delete $OWNER{$self};
    }
}

*swig_Visit_get = *Chemistry::OpenBabelc::OBMolAtomDFSIter_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBMolAtomDFSIter_Visit_set;
*Clear = *Chemistry::OpenBabelc::OBMolAtomDFSIter_Clear;
*SetIdx = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetIdx;
*SetHyb = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetHyb;
*SetAtomicNum = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetAtomicNum;
*SetIsotope = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetIsotope;
*SetImplicitValence = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetImplicitValence;
*IncrementImplicitValence = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IncrementImplicitValence;
*DecrementImplicitValence = *Chemistry::OpenBabelc::OBMolAtomDFSIter_DecrementImplicitValence;
*SetFormalCharge = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetFormalCharge;
*SetSpinMultiplicity = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetSpinMultiplicity;
*SetType = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetType;
*SetPartialCharge = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetPartialCharge;
*SetVector = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetVector;
*SetCoordPtr = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetCoordPtr;
*SetResidue = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetResidue;
*SetParent = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetParent;
*SetAromatic = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetAromatic;
*UnsetAromatic = *Chemistry::OpenBabelc::OBMolAtomDFSIter_UnsetAromatic;
*SetClockwiseStereo = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetClockwiseStereo;
*SetAntiClockwiseStereo = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetAntiClockwiseStereo;
*SetPositiveStereo = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetPositiveStereo;
*SetNegativeStereo = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetNegativeStereo;
*UnsetStereo = *Chemistry::OpenBabelc::OBMolAtomDFSIter_UnsetStereo;
*SetInRing = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetInRing;
*SetChiral = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetChiral;
*ClearCoordPtr = *Chemistry::OpenBabelc::OBMolAtomDFSIter_ClearCoordPtr;
*GetFormalCharge = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetFormalCharge;
*GetAtomicNum = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetAtomicNum;
*GetIsotope = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetIsotope;
*GetSpinMultiplicity = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetSpinMultiplicity;
*GetAtomicMass = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetAtomicMass;
*GetExactMass = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetExactMass;
*GetIdx = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetIdx;
*GetCoordinateIdx = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetCoordinateIdx;
*GetCIdx = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetCIdx;
*GetValence = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetValence;
*GetHyb = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetHyb;
*GetImplicitValence = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetImplicitValence;
*GetHvyValence = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetHvyValence;
*GetHeteroValence = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetHeteroValence;
*GetType = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetType;
*GetX = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetX;
*x = *Chemistry::OpenBabelc::OBMolAtomDFSIter_x;
*GetY = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetY;
*y = *Chemistry::OpenBabelc::OBMolAtomDFSIter_y;
*GetZ = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetZ;
*z = *Chemistry::OpenBabelc::OBMolAtomDFSIter_z;
*GetCoordinate = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetCoordinate;
*GetVector = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetVector;
*GetPartialCharge = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetPartialCharge;
*GetResidue = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetResidue;
*GetParent = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetParent;
*GetNewBondVector = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetNewBondVector;
*GetBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetBond;
*GetNextAtom = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetNextAtom;
*BeginBonds = *Chemistry::OpenBabelc::OBMolAtomDFSIter_BeginBonds;
*EndBonds = *Chemistry::OpenBabelc::OBMolAtomDFSIter_EndBonds;
*BeginBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_NextBond;
*BeginNbrAtom = *Chemistry::OpenBabelc::OBMolAtomDFSIter_BeginNbrAtom;
*NextNbrAtom = *Chemistry::OpenBabelc::OBMolAtomDFSIter_NextNbrAtom;
*GetDistance = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetDistance;
*GetAngle = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetAngle;
*NewResidue = *Chemistry::OpenBabelc::OBMolAtomDFSIter_NewResidue;
*DeleteResidue = *Chemistry::OpenBabelc::OBMolAtomDFSIter_DeleteResidue;
*AddBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_AddBond;
*InsertBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_InsertBond;
*DeleteBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_DeleteBond;
*ClearBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_ClearBond;
*CountFreeOxygens = *Chemistry::OpenBabelc::OBMolAtomDFSIter_CountFreeOxygens;
*ImplicitHydrogenCount = *Chemistry::OpenBabelc::OBMolAtomDFSIter_ImplicitHydrogenCount;
*ExplicitHydrogenCount = *Chemistry::OpenBabelc::OBMolAtomDFSIter_ExplicitHydrogenCount;
*MemberOfRingCount = *Chemistry::OpenBabelc::OBMolAtomDFSIter_MemberOfRingCount;
*MemberOfRingSize = *Chemistry::OpenBabelc::OBMolAtomDFSIter_MemberOfRingSize;
*CountRingBonds = *Chemistry::OpenBabelc::OBMolAtomDFSIter_CountRingBonds;
*SmallestBondAngle = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SmallestBondAngle;
*AverageBondAngle = *Chemistry::OpenBabelc::OBMolAtomDFSIter_AverageBondAngle;
*BOSum = *Chemistry::OpenBabelc::OBMolAtomDFSIter_BOSum;
*KBOSum = *Chemistry::OpenBabelc::OBMolAtomDFSIter_KBOSum;
*HtoMethyl = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HtoMethyl;
*SetHybAndGeom = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetHybAndGeom;
*ForceNoH = *Chemistry::OpenBabelc::OBMolAtomDFSIter_ForceNoH;
*HasNoHForced = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasNoHForced;
*HasResidue = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasResidue;
*IsHydrogen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsHydrogen;
*IsCarbon = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsCarbon;
*IsNitrogen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsNitrogen;
*IsOxygen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsOxygen;
*IsSulfur = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsSulfur;
*IsPhosphorus = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsPhosphorus;
*IsAromatic = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsInRing;
*IsInRingSize = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsInRingSize;
*IsHeteroatom = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsHeteroatom;
*IsNotCorH = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsNotCorH;
*IsConnected = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsConnected;
*IsOneThree = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsOneThree;
*IsOneFour = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsOneFour;
*IsCarboxylOxygen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsCarboxylOxygen;
*IsPhosphateOxygen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsPhosphateOxygen;
*IsSulfateOxygen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsSulfateOxygen;
*IsNitroOxygen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsNitroOxygen;
*IsAmideNitrogen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsAmideNitrogen;
*IsPolarHydrogen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsPolarHydrogen;
*IsNonPolarHydrogen = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsNonPolarHydrogen;
*IsAromaticNOxide = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsAromaticNOxide;
*IsChiral = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsChiral;
*IsAxial = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsAxial;
*IsClockwise = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsClockwise;
*IsAntiClockwise = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsAntiClockwise;
*IsPositiveStereo = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsPositiveStereo;
*IsNegativeStereo = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsNegativeStereo;
*HasChiralitySpecified = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasChiralitySpecified;
*HasChiralVolume = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasChiralVolume;
*IsHbondAcceptor = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsHbondAcceptor;
*IsHbondDonor = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsHbondDonor;
*IsHbondDonorH = *Chemistry::OpenBabelc::OBMolAtomDFSIter_IsHbondDonorH;
*HasAlphaBetaUnsat = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasAlphaBetaUnsat;
*HasBondOfOrder = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasBondOfOrder;
*CountBondsOfOrder = *Chemistry::OpenBabelc::OBMolAtomDFSIter_CountBondsOfOrder;
*HasNonSingleBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasNonSingleBond;
*HasSingleBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasSingleBond;
*HasDoubleBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasDoubleBond;
*HasAromaticBond = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasAromaticBond;
*MatchesSMARTS = *Chemistry::OpenBabelc::OBMolAtomDFSIter_MatchesSMARTS;
*DoTransformations = *Chemistry::OpenBabelc::OBMolAtomDFSIter_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBMolAtomDFSIter_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBMolAtomDFSIter_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBMolAtomDFSIter_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBMolAtomDFSIter_SetData;
*DataSize = *Chemistry::OpenBabelc::OBMolAtomDFSIter_DataSize;
*GetData = *Chemistry::OpenBabelc::OBMolAtomDFSIter_GetData;
*BeginData = *Chemistry::OpenBabelc::OBMolAtomDFSIter_BeginData;
*EndData = *Chemistry::OpenBabelc::OBMolAtomDFSIter_EndData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBMolAtomBFSIter ##############

package Chemistry::OpenBabel::OBMolAtomBFSIter;
use overload
    "++" => sub { $_[0]->__plusplus__()},
    "fallback" => 1;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBMolAtomBFSIter(@_);
    bless $self, $pkg if defined($self);
}

*__plusplus__ = *Chemistry::OpenBabelc::OBMolAtomBFSIter___plusplus__;
*__deref__ = *Chemistry::OpenBabelc::OBMolAtomBFSIter___deref__;
*__ref__ = *Chemistry::OpenBabelc::OBMolAtomBFSIter___ref__;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBMolAtomBFSIter($self);
        delete $OWNER{$self};
    }
}

*swig_Visit_get = *Chemistry::OpenBabelc::OBMolAtomBFSIter_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBMolAtomBFSIter_Visit_set;
*Clear = *Chemistry::OpenBabelc::OBMolAtomBFSIter_Clear;
*SetIdx = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetIdx;
*SetHyb = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetHyb;
*SetAtomicNum = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetAtomicNum;
*SetIsotope = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetIsotope;
*SetImplicitValence = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetImplicitValence;
*IncrementImplicitValence = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IncrementImplicitValence;
*DecrementImplicitValence = *Chemistry::OpenBabelc::OBMolAtomBFSIter_DecrementImplicitValence;
*SetFormalCharge = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetFormalCharge;
*SetSpinMultiplicity = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetSpinMultiplicity;
*SetType = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetType;
*SetPartialCharge = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetPartialCharge;
*SetVector = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetVector;
*SetCoordPtr = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetCoordPtr;
*SetResidue = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetResidue;
*SetParent = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetParent;
*SetAromatic = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetAromatic;
*UnsetAromatic = *Chemistry::OpenBabelc::OBMolAtomBFSIter_UnsetAromatic;
*SetClockwiseStereo = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetClockwiseStereo;
*SetAntiClockwiseStereo = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetAntiClockwiseStereo;
*SetPositiveStereo = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetPositiveStereo;
*SetNegativeStereo = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetNegativeStereo;
*UnsetStereo = *Chemistry::OpenBabelc::OBMolAtomBFSIter_UnsetStereo;
*SetInRing = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetInRing;
*SetChiral = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetChiral;
*ClearCoordPtr = *Chemistry::OpenBabelc::OBMolAtomBFSIter_ClearCoordPtr;
*GetFormalCharge = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetFormalCharge;
*GetAtomicNum = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetAtomicNum;
*GetIsotope = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetIsotope;
*GetSpinMultiplicity = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetSpinMultiplicity;
*GetAtomicMass = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetAtomicMass;
*GetExactMass = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetExactMass;
*GetIdx = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetIdx;
*GetCoordinateIdx = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetCoordinateIdx;
*GetCIdx = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetCIdx;
*GetValence = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetValence;
*GetHyb = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetHyb;
*GetImplicitValence = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetImplicitValence;
*GetHvyValence = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetHvyValence;
*GetHeteroValence = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetHeteroValence;
*GetType = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetType;
*GetX = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetX;
*x = *Chemistry::OpenBabelc::OBMolAtomBFSIter_x;
*GetY = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetY;
*y = *Chemistry::OpenBabelc::OBMolAtomBFSIter_y;
*GetZ = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetZ;
*z = *Chemistry::OpenBabelc::OBMolAtomBFSIter_z;
*GetCoordinate = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetCoordinate;
*GetVector = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetVector;
*GetPartialCharge = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetPartialCharge;
*GetResidue = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetResidue;
*GetParent = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetParent;
*GetNewBondVector = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetNewBondVector;
*GetBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetBond;
*GetNextAtom = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetNextAtom;
*BeginBonds = *Chemistry::OpenBabelc::OBMolAtomBFSIter_BeginBonds;
*EndBonds = *Chemistry::OpenBabelc::OBMolAtomBFSIter_EndBonds;
*BeginBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_NextBond;
*BeginNbrAtom = *Chemistry::OpenBabelc::OBMolAtomBFSIter_BeginNbrAtom;
*NextNbrAtom = *Chemistry::OpenBabelc::OBMolAtomBFSIter_NextNbrAtom;
*GetDistance = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetDistance;
*GetAngle = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetAngle;
*NewResidue = *Chemistry::OpenBabelc::OBMolAtomBFSIter_NewResidue;
*DeleteResidue = *Chemistry::OpenBabelc::OBMolAtomBFSIter_DeleteResidue;
*AddBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_AddBond;
*InsertBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_InsertBond;
*DeleteBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_DeleteBond;
*ClearBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_ClearBond;
*CountFreeOxygens = *Chemistry::OpenBabelc::OBMolAtomBFSIter_CountFreeOxygens;
*ImplicitHydrogenCount = *Chemistry::OpenBabelc::OBMolAtomBFSIter_ImplicitHydrogenCount;
*ExplicitHydrogenCount = *Chemistry::OpenBabelc::OBMolAtomBFSIter_ExplicitHydrogenCount;
*MemberOfRingCount = *Chemistry::OpenBabelc::OBMolAtomBFSIter_MemberOfRingCount;
*MemberOfRingSize = *Chemistry::OpenBabelc::OBMolAtomBFSIter_MemberOfRingSize;
*CountRingBonds = *Chemistry::OpenBabelc::OBMolAtomBFSIter_CountRingBonds;
*SmallestBondAngle = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SmallestBondAngle;
*AverageBondAngle = *Chemistry::OpenBabelc::OBMolAtomBFSIter_AverageBondAngle;
*BOSum = *Chemistry::OpenBabelc::OBMolAtomBFSIter_BOSum;
*KBOSum = *Chemistry::OpenBabelc::OBMolAtomBFSIter_KBOSum;
*HtoMethyl = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HtoMethyl;
*SetHybAndGeom = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetHybAndGeom;
*ForceNoH = *Chemistry::OpenBabelc::OBMolAtomBFSIter_ForceNoH;
*HasNoHForced = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasNoHForced;
*HasResidue = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasResidue;
*IsHydrogen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsHydrogen;
*IsCarbon = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsCarbon;
*IsNitrogen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsNitrogen;
*IsOxygen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsOxygen;
*IsSulfur = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsSulfur;
*IsPhosphorus = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsPhosphorus;
*IsAromatic = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsInRing;
*IsInRingSize = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsInRingSize;
*IsHeteroatom = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsHeteroatom;
*IsNotCorH = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsNotCorH;
*IsConnected = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsConnected;
*IsOneThree = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsOneThree;
*IsOneFour = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsOneFour;
*IsCarboxylOxygen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsCarboxylOxygen;
*IsPhosphateOxygen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsPhosphateOxygen;
*IsSulfateOxygen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsSulfateOxygen;
*IsNitroOxygen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsNitroOxygen;
*IsAmideNitrogen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsAmideNitrogen;
*IsPolarHydrogen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsPolarHydrogen;
*IsNonPolarHydrogen = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsNonPolarHydrogen;
*IsAromaticNOxide = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsAromaticNOxide;
*IsChiral = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsChiral;
*IsAxial = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsAxial;
*IsClockwise = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsClockwise;
*IsAntiClockwise = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsAntiClockwise;
*IsPositiveStereo = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsPositiveStereo;
*IsNegativeStereo = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsNegativeStereo;
*HasChiralitySpecified = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasChiralitySpecified;
*HasChiralVolume = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasChiralVolume;
*IsHbondAcceptor = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsHbondAcceptor;
*IsHbondDonor = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsHbondDonor;
*IsHbondDonorH = *Chemistry::OpenBabelc::OBMolAtomBFSIter_IsHbondDonorH;
*HasAlphaBetaUnsat = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasAlphaBetaUnsat;
*HasBondOfOrder = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasBondOfOrder;
*CountBondsOfOrder = *Chemistry::OpenBabelc::OBMolAtomBFSIter_CountBondsOfOrder;
*HasNonSingleBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasNonSingleBond;
*HasSingleBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasSingleBond;
*HasDoubleBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasDoubleBond;
*HasAromaticBond = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasAromaticBond;
*MatchesSMARTS = *Chemistry::OpenBabelc::OBMolAtomBFSIter_MatchesSMARTS;
*DoTransformations = *Chemistry::OpenBabelc::OBMolAtomBFSIter_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBMolAtomBFSIter_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBMolAtomBFSIter_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBMolAtomBFSIter_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBMolAtomBFSIter_SetData;
*DataSize = *Chemistry::OpenBabelc::OBMolAtomBFSIter_DataSize;
*GetData = *Chemistry::OpenBabelc::OBMolAtomBFSIter_GetData;
*BeginData = *Chemistry::OpenBabelc::OBMolAtomBFSIter_BeginData;
*EndData = *Chemistry::OpenBabelc::OBMolAtomBFSIter_EndData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBMolBondIter ##############

package Chemistry::OpenBabel::OBMolBondIter;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBMolBondIter(@_);
    bless $self, $pkg if defined($self);
}

*good = *Chemistry::OpenBabelc::OBMolBondIter_good;
*inc = *Chemistry::OpenBabelc::OBMolBondIter_inc;
*deref = *Chemistry::OpenBabelc::OBMolBondIter_deref;
*__ref__ = *Chemistry::OpenBabelc::OBMolBondIter___ref__;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBMolBondIter($self);
        delete $OWNER{$self};
    }
}

*swig_Visit_get = *Chemistry::OpenBabelc::OBMolBondIter_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBMolBondIter_Visit_set;
*SetIdx = *Chemistry::OpenBabelc::OBMolBondIter_SetIdx;
*SetBO = *Chemistry::OpenBabelc::OBMolBondIter_SetBO;
*SetBegin = *Chemistry::OpenBabelc::OBMolBondIter_SetBegin;
*SetEnd = *Chemistry::OpenBabelc::OBMolBondIter_SetEnd;
*SetParent = *Chemistry::OpenBabelc::OBMolBondIter_SetParent;
*SetLength = *Chemistry::OpenBabelc::OBMolBondIter_SetLength;
*Set = *Chemistry::OpenBabelc::OBMolBondIter_Set;
*SetKSingle = *Chemistry::OpenBabelc::OBMolBondIter_SetKSingle;
*SetKDouble = *Chemistry::OpenBabelc::OBMolBondIter_SetKDouble;
*SetKTriple = *Chemistry::OpenBabelc::OBMolBondIter_SetKTriple;
*SetAromatic = *Chemistry::OpenBabelc::OBMolBondIter_SetAromatic;
*SetHash = *Chemistry::OpenBabelc::OBMolBondIter_SetHash;
*SetWedge = *Chemistry::OpenBabelc::OBMolBondIter_SetWedge;
*SetUp = *Chemistry::OpenBabelc::OBMolBondIter_SetUp;
*SetDown = *Chemistry::OpenBabelc::OBMolBondIter_SetDown;
*SetInRing = *Chemistry::OpenBabelc::OBMolBondIter_SetInRing;
*SetClosure = *Chemistry::OpenBabelc::OBMolBondIter_SetClosure;
*UnsetHash = *Chemistry::OpenBabelc::OBMolBondIter_UnsetHash;
*UnsetWedge = *Chemistry::OpenBabelc::OBMolBondIter_UnsetWedge;
*UnsetUp = *Chemistry::OpenBabelc::OBMolBondIter_UnsetUp;
*UnsetDown = *Chemistry::OpenBabelc::OBMolBondIter_UnsetDown;
*UnsetAromatic = *Chemistry::OpenBabelc::OBMolBondIter_UnsetAromatic;
*UnsetKekule = *Chemistry::OpenBabelc::OBMolBondIter_UnsetKekule;
*GetIdx = *Chemistry::OpenBabelc::OBMolBondIter_GetIdx;
*GetBO = *Chemistry::OpenBabelc::OBMolBondIter_GetBO;
*GetBondOrder = *Chemistry::OpenBabelc::OBMolBondIter_GetBondOrder;
*GetFlags = *Chemistry::OpenBabelc::OBMolBondIter_GetFlags;
*GetBeginAtomIdx = *Chemistry::OpenBabelc::OBMolBondIter_GetBeginAtomIdx;
*GetEndAtomIdx = *Chemistry::OpenBabelc::OBMolBondIter_GetEndAtomIdx;
*GetBeginAtom = *Chemistry::OpenBabelc::OBMolBondIter_GetBeginAtom;
*GetEndAtom = *Chemistry::OpenBabelc::OBMolBondIter_GetEndAtom;
*GetNbrAtom = *Chemistry::OpenBabelc::OBMolBondIter_GetNbrAtom;
*GetParent = *Chemistry::OpenBabelc::OBMolBondIter_GetParent;
*GetEquibLength = *Chemistry::OpenBabelc::OBMolBondIter_GetEquibLength;
*GetLength = *Chemistry::OpenBabelc::OBMolBondIter_GetLength;
*GetNbrAtomIdx = *Chemistry::OpenBabelc::OBMolBondIter_GetNbrAtomIdx;
*IsAromatic = *Chemistry::OpenBabelc::OBMolBondIter_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBMolBondIter_IsInRing;
*IsRotor = *Chemistry::OpenBabelc::OBMolBondIter_IsRotor;
*IsAmide = *Chemistry::OpenBabelc::OBMolBondIter_IsAmide;
*IsPrimaryAmide = *Chemistry::OpenBabelc::OBMolBondIter_IsPrimaryAmide;
*IsSecondaryAmide = *Chemistry::OpenBabelc::OBMolBondIter_IsSecondaryAmide;
*IsEster = *Chemistry::OpenBabelc::OBMolBondIter_IsEster;
*IsCarbonyl = *Chemistry::OpenBabelc::OBMolBondIter_IsCarbonyl;
*IsSingle = *Chemistry::OpenBabelc::OBMolBondIter_IsSingle;
*IsDouble = *Chemistry::OpenBabelc::OBMolBondIter_IsDouble;
*IsTriple = *Chemistry::OpenBabelc::OBMolBondIter_IsTriple;
*IsKSingle = *Chemistry::OpenBabelc::OBMolBondIter_IsKSingle;
*IsKDouble = *Chemistry::OpenBabelc::OBMolBondIter_IsKDouble;
*IsKTriple = *Chemistry::OpenBabelc::OBMolBondIter_IsKTriple;
*IsClosure = *Chemistry::OpenBabelc::OBMolBondIter_IsClosure;
*IsUp = *Chemistry::OpenBabelc::OBMolBondIter_IsUp;
*IsDown = *Chemistry::OpenBabelc::OBMolBondIter_IsDown;
*IsWedge = *Chemistry::OpenBabelc::OBMolBondIter_IsWedge;
*IsHash = *Chemistry::OpenBabelc::OBMolBondIter_IsHash;
*IsDoubleBondGeometry = *Chemistry::OpenBabelc::OBMolBondIter_IsDoubleBondGeometry;
*DoTransformations = *Chemistry::OpenBabelc::OBMolBondIter_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBMolBondIter_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBMolBondIter_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBMolBondIter_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBMolBondIter_SetData;
*DataSize = *Chemistry::OpenBabelc::OBMolBondIter_DataSize;
*GetData = *Chemistry::OpenBabelc::OBMolBondIter_GetData;
*BeginData = *Chemistry::OpenBabelc::OBMolBondIter_BeginData;
*EndData = *Chemistry::OpenBabelc::OBMolBondIter_EndData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBAtomAtomIter ##############

package Chemistry::OpenBabel::OBAtomAtomIter;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBAtomAtomIter(@_);
    bless $self, $pkg if defined($self);
}

*good = *Chemistry::OpenBabelc::OBAtomAtomIter_good;
*inc = *Chemistry::OpenBabelc::OBAtomAtomIter_inc;
*deref = *Chemistry::OpenBabelc::OBAtomAtomIter_deref;
*__ref__ = *Chemistry::OpenBabelc::OBAtomAtomIter___ref__;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBAtomAtomIter($self);
        delete $OWNER{$self};
    }
}

*swig_Visit_get = *Chemistry::OpenBabelc::OBAtomAtomIter_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBAtomAtomIter_Visit_set;
*Clear = *Chemistry::OpenBabelc::OBAtomAtomIter_Clear;
*SetIdx = *Chemistry::OpenBabelc::OBAtomAtomIter_SetIdx;
*SetHyb = *Chemistry::OpenBabelc::OBAtomAtomIter_SetHyb;
*SetAtomicNum = *Chemistry::OpenBabelc::OBAtomAtomIter_SetAtomicNum;
*SetIsotope = *Chemistry::OpenBabelc::OBAtomAtomIter_SetIsotope;
*SetImplicitValence = *Chemistry::OpenBabelc::OBAtomAtomIter_SetImplicitValence;
*IncrementImplicitValence = *Chemistry::OpenBabelc::OBAtomAtomIter_IncrementImplicitValence;
*DecrementImplicitValence = *Chemistry::OpenBabelc::OBAtomAtomIter_DecrementImplicitValence;
*SetFormalCharge = *Chemistry::OpenBabelc::OBAtomAtomIter_SetFormalCharge;
*SetSpinMultiplicity = *Chemistry::OpenBabelc::OBAtomAtomIter_SetSpinMultiplicity;
*SetType = *Chemistry::OpenBabelc::OBAtomAtomIter_SetType;
*SetPartialCharge = *Chemistry::OpenBabelc::OBAtomAtomIter_SetPartialCharge;
*SetVector = *Chemistry::OpenBabelc::OBAtomAtomIter_SetVector;
*SetCoordPtr = *Chemistry::OpenBabelc::OBAtomAtomIter_SetCoordPtr;
*SetResidue = *Chemistry::OpenBabelc::OBAtomAtomIter_SetResidue;
*SetParent = *Chemistry::OpenBabelc::OBAtomAtomIter_SetParent;
*SetAromatic = *Chemistry::OpenBabelc::OBAtomAtomIter_SetAromatic;
*UnsetAromatic = *Chemistry::OpenBabelc::OBAtomAtomIter_UnsetAromatic;
*SetClockwiseStereo = *Chemistry::OpenBabelc::OBAtomAtomIter_SetClockwiseStereo;
*SetAntiClockwiseStereo = *Chemistry::OpenBabelc::OBAtomAtomIter_SetAntiClockwiseStereo;
*SetPositiveStereo = *Chemistry::OpenBabelc::OBAtomAtomIter_SetPositiveStereo;
*SetNegativeStereo = *Chemistry::OpenBabelc::OBAtomAtomIter_SetNegativeStereo;
*UnsetStereo = *Chemistry::OpenBabelc::OBAtomAtomIter_UnsetStereo;
*SetInRing = *Chemistry::OpenBabelc::OBAtomAtomIter_SetInRing;
*SetChiral = *Chemistry::OpenBabelc::OBAtomAtomIter_SetChiral;
*ClearCoordPtr = *Chemistry::OpenBabelc::OBAtomAtomIter_ClearCoordPtr;
*GetFormalCharge = *Chemistry::OpenBabelc::OBAtomAtomIter_GetFormalCharge;
*GetAtomicNum = *Chemistry::OpenBabelc::OBAtomAtomIter_GetAtomicNum;
*GetIsotope = *Chemistry::OpenBabelc::OBAtomAtomIter_GetIsotope;
*GetSpinMultiplicity = *Chemistry::OpenBabelc::OBAtomAtomIter_GetSpinMultiplicity;
*GetAtomicMass = *Chemistry::OpenBabelc::OBAtomAtomIter_GetAtomicMass;
*GetExactMass = *Chemistry::OpenBabelc::OBAtomAtomIter_GetExactMass;
*GetIdx = *Chemistry::OpenBabelc::OBAtomAtomIter_GetIdx;
*GetCoordinateIdx = *Chemistry::OpenBabelc::OBAtomAtomIter_GetCoordinateIdx;
*GetCIdx = *Chemistry::OpenBabelc::OBAtomAtomIter_GetCIdx;
*GetValence = *Chemistry::OpenBabelc::OBAtomAtomIter_GetValence;
*GetHyb = *Chemistry::OpenBabelc::OBAtomAtomIter_GetHyb;
*GetImplicitValence = *Chemistry::OpenBabelc::OBAtomAtomIter_GetImplicitValence;
*GetHvyValence = *Chemistry::OpenBabelc::OBAtomAtomIter_GetHvyValence;
*GetHeteroValence = *Chemistry::OpenBabelc::OBAtomAtomIter_GetHeteroValence;
*GetType = *Chemistry::OpenBabelc::OBAtomAtomIter_GetType;
*GetX = *Chemistry::OpenBabelc::OBAtomAtomIter_GetX;
*x = *Chemistry::OpenBabelc::OBAtomAtomIter_x;
*GetY = *Chemistry::OpenBabelc::OBAtomAtomIter_GetY;
*y = *Chemistry::OpenBabelc::OBAtomAtomIter_y;
*GetZ = *Chemistry::OpenBabelc::OBAtomAtomIter_GetZ;
*z = *Chemistry::OpenBabelc::OBAtomAtomIter_z;
*GetCoordinate = *Chemistry::OpenBabelc::OBAtomAtomIter_GetCoordinate;
*GetVector = *Chemistry::OpenBabelc::OBAtomAtomIter_GetVector;
*GetPartialCharge = *Chemistry::OpenBabelc::OBAtomAtomIter_GetPartialCharge;
*GetResidue = *Chemistry::OpenBabelc::OBAtomAtomIter_GetResidue;
*GetParent = *Chemistry::OpenBabelc::OBAtomAtomIter_GetParent;
*GetNewBondVector = *Chemistry::OpenBabelc::OBAtomAtomIter_GetNewBondVector;
*GetBond = *Chemistry::OpenBabelc::OBAtomAtomIter_GetBond;
*GetNextAtom = *Chemistry::OpenBabelc::OBAtomAtomIter_GetNextAtom;
*BeginBonds = *Chemistry::OpenBabelc::OBAtomAtomIter_BeginBonds;
*EndBonds = *Chemistry::OpenBabelc::OBAtomAtomIter_EndBonds;
*BeginBond = *Chemistry::OpenBabelc::OBAtomAtomIter_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBAtomAtomIter_NextBond;
*BeginNbrAtom = *Chemistry::OpenBabelc::OBAtomAtomIter_BeginNbrAtom;
*NextNbrAtom = *Chemistry::OpenBabelc::OBAtomAtomIter_NextNbrAtom;
*GetDistance = *Chemistry::OpenBabelc::OBAtomAtomIter_GetDistance;
*GetAngle = *Chemistry::OpenBabelc::OBAtomAtomIter_GetAngle;
*NewResidue = *Chemistry::OpenBabelc::OBAtomAtomIter_NewResidue;
*DeleteResidue = *Chemistry::OpenBabelc::OBAtomAtomIter_DeleteResidue;
*AddBond = *Chemistry::OpenBabelc::OBAtomAtomIter_AddBond;
*InsertBond = *Chemistry::OpenBabelc::OBAtomAtomIter_InsertBond;
*DeleteBond = *Chemistry::OpenBabelc::OBAtomAtomIter_DeleteBond;
*ClearBond = *Chemistry::OpenBabelc::OBAtomAtomIter_ClearBond;
*CountFreeOxygens = *Chemistry::OpenBabelc::OBAtomAtomIter_CountFreeOxygens;
*ImplicitHydrogenCount = *Chemistry::OpenBabelc::OBAtomAtomIter_ImplicitHydrogenCount;
*ExplicitHydrogenCount = *Chemistry::OpenBabelc::OBAtomAtomIter_ExplicitHydrogenCount;
*MemberOfRingCount = *Chemistry::OpenBabelc::OBAtomAtomIter_MemberOfRingCount;
*MemberOfRingSize = *Chemistry::OpenBabelc::OBAtomAtomIter_MemberOfRingSize;
*CountRingBonds = *Chemistry::OpenBabelc::OBAtomAtomIter_CountRingBonds;
*SmallestBondAngle = *Chemistry::OpenBabelc::OBAtomAtomIter_SmallestBondAngle;
*AverageBondAngle = *Chemistry::OpenBabelc::OBAtomAtomIter_AverageBondAngle;
*BOSum = *Chemistry::OpenBabelc::OBAtomAtomIter_BOSum;
*KBOSum = *Chemistry::OpenBabelc::OBAtomAtomIter_KBOSum;
*HtoMethyl = *Chemistry::OpenBabelc::OBAtomAtomIter_HtoMethyl;
*SetHybAndGeom = *Chemistry::OpenBabelc::OBAtomAtomIter_SetHybAndGeom;
*ForceNoH = *Chemistry::OpenBabelc::OBAtomAtomIter_ForceNoH;
*HasNoHForced = *Chemistry::OpenBabelc::OBAtomAtomIter_HasNoHForced;
*HasResidue = *Chemistry::OpenBabelc::OBAtomAtomIter_HasResidue;
*IsHydrogen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsHydrogen;
*IsCarbon = *Chemistry::OpenBabelc::OBAtomAtomIter_IsCarbon;
*IsNitrogen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsNitrogen;
*IsOxygen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsOxygen;
*IsSulfur = *Chemistry::OpenBabelc::OBAtomAtomIter_IsSulfur;
*IsPhosphorus = *Chemistry::OpenBabelc::OBAtomAtomIter_IsPhosphorus;
*IsAromatic = *Chemistry::OpenBabelc::OBAtomAtomIter_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBAtomAtomIter_IsInRing;
*IsInRingSize = *Chemistry::OpenBabelc::OBAtomAtomIter_IsInRingSize;
*IsHeteroatom = *Chemistry::OpenBabelc::OBAtomAtomIter_IsHeteroatom;
*IsNotCorH = *Chemistry::OpenBabelc::OBAtomAtomIter_IsNotCorH;
*IsConnected = *Chemistry::OpenBabelc::OBAtomAtomIter_IsConnected;
*IsOneThree = *Chemistry::OpenBabelc::OBAtomAtomIter_IsOneThree;
*IsOneFour = *Chemistry::OpenBabelc::OBAtomAtomIter_IsOneFour;
*IsCarboxylOxygen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsCarboxylOxygen;
*IsPhosphateOxygen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsPhosphateOxygen;
*IsSulfateOxygen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsSulfateOxygen;
*IsNitroOxygen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsNitroOxygen;
*IsAmideNitrogen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsAmideNitrogen;
*IsPolarHydrogen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsPolarHydrogen;
*IsNonPolarHydrogen = *Chemistry::OpenBabelc::OBAtomAtomIter_IsNonPolarHydrogen;
*IsAromaticNOxide = *Chemistry::OpenBabelc::OBAtomAtomIter_IsAromaticNOxide;
*IsChiral = *Chemistry::OpenBabelc::OBAtomAtomIter_IsChiral;
*IsAxial = *Chemistry::OpenBabelc::OBAtomAtomIter_IsAxial;
*IsClockwise = *Chemistry::OpenBabelc::OBAtomAtomIter_IsClockwise;
*IsAntiClockwise = *Chemistry::OpenBabelc::OBAtomAtomIter_IsAntiClockwise;
*IsPositiveStereo = *Chemistry::OpenBabelc::OBAtomAtomIter_IsPositiveStereo;
*IsNegativeStereo = *Chemistry::OpenBabelc::OBAtomAtomIter_IsNegativeStereo;
*HasChiralitySpecified = *Chemistry::OpenBabelc::OBAtomAtomIter_HasChiralitySpecified;
*HasChiralVolume = *Chemistry::OpenBabelc::OBAtomAtomIter_HasChiralVolume;
*IsHbondAcceptor = *Chemistry::OpenBabelc::OBAtomAtomIter_IsHbondAcceptor;
*IsHbondDonor = *Chemistry::OpenBabelc::OBAtomAtomIter_IsHbondDonor;
*IsHbondDonorH = *Chemistry::OpenBabelc::OBAtomAtomIter_IsHbondDonorH;
*HasAlphaBetaUnsat = *Chemistry::OpenBabelc::OBAtomAtomIter_HasAlphaBetaUnsat;
*HasBondOfOrder = *Chemistry::OpenBabelc::OBAtomAtomIter_HasBondOfOrder;
*CountBondsOfOrder = *Chemistry::OpenBabelc::OBAtomAtomIter_CountBondsOfOrder;
*HasNonSingleBond = *Chemistry::OpenBabelc::OBAtomAtomIter_HasNonSingleBond;
*HasSingleBond = *Chemistry::OpenBabelc::OBAtomAtomIter_HasSingleBond;
*HasDoubleBond = *Chemistry::OpenBabelc::OBAtomAtomIter_HasDoubleBond;
*HasAromaticBond = *Chemistry::OpenBabelc::OBAtomAtomIter_HasAromaticBond;
*MatchesSMARTS = *Chemistry::OpenBabelc::OBAtomAtomIter_MatchesSMARTS;
*DoTransformations = *Chemistry::OpenBabelc::OBAtomAtomIter_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBAtomAtomIter_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBAtomAtomIter_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBAtomAtomIter_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBAtomAtomIter_SetData;
*DataSize = *Chemistry::OpenBabelc::OBAtomAtomIter_DataSize;
*GetData = *Chemistry::OpenBabelc::OBAtomAtomIter_GetData;
*BeginData = *Chemistry::OpenBabelc::OBAtomAtomIter_BeginData;
*EndData = *Chemistry::OpenBabelc::OBAtomAtomIter_EndData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBAtomBondIter ##############

package Chemistry::OpenBabel::OBAtomBondIter;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBAtomBondIter(@_);
    bless $self, $pkg if defined($self);
}

*good = *Chemistry::OpenBabelc::OBAtomBondIter_good;
*inc = *Chemistry::OpenBabelc::OBAtomBondIter_inc;
*deref = *Chemistry::OpenBabelc::OBAtomBondIter_deref;
*__ref__ = *Chemistry::OpenBabelc::OBAtomBondIter___ref__;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBAtomBondIter($self);
        delete $OWNER{$self};
    }
}

*swig_Visit_get = *Chemistry::OpenBabelc::OBAtomBondIter_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBAtomBondIter_Visit_set;
*SetIdx = *Chemistry::OpenBabelc::OBAtomBondIter_SetIdx;
*SetBO = *Chemistry::OpenBabelc::OBAtomBondIter_SetBO;
*SetBegin = *Chemistry::OpenBabelc::OBAtomBondIter_SetBegin;
*SetEnd = *Chemistry::OpenBabelc::OBAtomBondIter_SetEnd;
*SetParent = *Chemistry::OpenBabelc::OBAtomBondIter_SetParent;
*SetLength = *Chemistry::OpenBabelc::OBAtomBondIter_SetLength;
*Set = *Chemistry::OpenBabelc::OBAtomBondIter_Set;
*SetKSingle = *Chemistry::OpenBabelc::OBAtomBondIter_SetKSingle;
*SetKDouble = *Chemistry::OpenBabelc::OBAtomBondIter_SetKDouble;
*SetKTriple = *Chemistry::OpenBabelc::OBAtomBondIter_SetKTriple;
*SetAromatic = *Chemistry::OpenBabelc::OBAtomBondIter_SetAromatic;
*SetHash = *Chemistry::OpenBabelc::OBAtomBondIter_SetHash;
*SetWedge = *Chemistry::OpenBabelc::OBAtomBondIter_SetWedge;
*SetUp = *Chemistry::OpenBabelc::OBAtomBondIter_SetUp;
*SetDown = *Chemistry::OpenBabelc::OBAtomBondIter_SetDown;
*SetInRing = *Chemistry::OpenBabelc::OBAtomBondIter_SetInRing;
*SetClosure = *Chemistry::OpenBabelc::OBAtomBondIter_SetClosure;
*UnsetHash = *Chemistry::OpenBabelc::OBAtomBondIter_UnsetHash;
*UnsetWedge = *Chemistry::OpenBabelc::OBAtomBondIter_UnsetWedge;
*UnsetUp = *Chemistry::OpenBabelc::OBAtomBondIter_UnsetUp;
*UnsetDown = *Chemistry::OpenBabelc::OBAtomBondIter_UnsetDown;
*UnsetAromatic = *Chemistry::OpenBabelc::OBAtomBondIter_UnsetAromatic;
*UnsetKekule = *Chemistry::OpenBabelc::OBAtomBondIter_UnsetKekule;
*GetIdx = *Chemistry::OpenBabelc::OBAtomBondIter_GetIdx;
*GetBO = *Chemistry::OpenBabelc::OBAtomBondIter_GetBO;
*GetBondOrder = *Chemistry::OpenBabelc::OBAtomBondIter_GetBondOrder;
*GetFlags = *Chemistry::OpenBabelc::OBAtomBondIter_GetFlags;
*GetBeginAtomIdx = *Chemistry::OpenBabelc::OBAtomBondIter_GetBeginAtomIdx;
*GetEndAtomIdx = *Chemistry::OpenBabelc::OBAtomBondIter_GetEndAtomIdx;
*GetBeginAtom = *Chemistry::OpenBabelc::OBAtomBondIter_GetBeginAtom;
*GetEndAtom = *Chemistry::OpenBabelc::OBAtomBondIter_GetEndAtom;
*GetNbrAtom = *Chemistry::OpenBabelc::OBAtomBondIter_GetNbrAtom;
*GetParent = *Chemistry::OpenBabelc::OBAtomBondIter_GetParent;
*GetEquibLength = *Chemistry::OpenBabelc::OBAtomBondIter_GetEquibLength;
*GetLength = *Chemistry::OpenBabelc::OBAtomBondIter_GetLength;
*GetNbrAtomIdx = *Chemistry::OpenBabelc::OBAtomBondIter_GetNbrAtomIdx;
*IsAromatic = *Chemistry::OpenBabelc::OBAtomBondIter_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBAtomBondIter_IsInRing;
*IsRotor = *Chemistry::OpenBabelc::OBAtomBondIter_IsRotor;
*IsAmide = *Chemistry::OpenBabelc::OBAtomBondIter_IsAmide;
*IsPrimaryAmide = *Chemistry::OpenBabelc::OBAtomBondIter_IsPrimaryAmide;
*IsSecondaryAmide = *Chemistry::OpenBabelc::OBAtomBondIter_IsSecondaryAmide;
*IsEster = *Chemistry::OpenBabelc::OBAtomBondIter_IsEster;
*IsCarbonyl = *Chemistry::OpenBabelc::OBAtomBondIter_IsCarbonyl;
*IsSingle = *Chemistry::OpenBabelc::OBAtomBondIter_IsSingle;
*IsDouble = *Chemistry::OpenBabelc::OBAtomBondIter_IsDouble;
*IsTriple = *Chemistry::OpenBabelc::OBAtomBondIter_IsTriple;
*IsKSingle = *Chemistry::OpenBabelc::OBAtomBondIter_IsKSingle;
*IsKDouble = *Chemistry::OpenBabelc::OBAtomBondIter_IsKDouble;
*IsKTriple = *Chemistry::OpenBabelc::OBAtomBondIter_IsKTriple;
*IsClosure = *Chemistry::OpenBabelc::OBAtomBondIter_IsClosure;
*IsUp = *Chemistry::OpenBabelc::OBAtomBondIter_IsUp;
*IsDown = *Chemistry::OpenBabelc::OBAtomBondIter_IsDown;
*IsWedge = *Chemistry::OpenBabelc::OBAtomBondIter_IsWedge;
*IsHash = *Chemistry::OpenBabelc::OBAtomBondIter_IsHash;
*IsDoubleBondGeometry = *Chemistry::OpenBabelc::OBAtomBondIter_IsDoubleBondGeometry;
*DoTransformations = *Chemistry::OpenBabelc::OBAtomBondIter_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBAtomBondIter_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBAtomBondIter_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBAtomBondIter_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBAtomBondIter_SetData;
*DataSize = *Chemistry::OpenBabelc::OBAtomBondIter_DataSize;
*GetData = *Chemistry::OpenBabelc::OBAtomBondIter_GetData;
*BeginData = *Chemistry::OpenBabelc::OBAtomBondIter_BeginData;
*EndData = *Chemistry::OpenBabelc::OBAtomBondIter_EndData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBResidueIter ##############

package Chemistry::OpenBabel::OBResidueIter;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBResidueIter(@_);
    bless $self, $pkg if defined($self);
}

*good = *Chemistry::OpenBabelc::OBResidueIter_good;
*inc = *Chemistry::OpenBabelc::OBResidueIter_inc;
*__deref__ = *Chemistry::OpenBabelc::OBResidueIter___deref__;
*__ref__ = *Chemistry::OpenBabelc::OBResidueIter___ref__;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBResidueIter($self);
        delete $OWNER{$self};
    }
}

*AddAtom = *Chemistry::OpenBabelc::OBResidueIter_AddAtom;
*InsertAtom = *Chemistry::OpenBabelc::OBResidueIter_InsertAtom;
*RemoveAtom = *Chemistry::OpenBabelc::OBResidueIter_RemoveAtom;
*Clear = *Chemistry::OpenBabelc::OBResidueIter_Clear;
*SetName = *Chemistry::OpenBabelc::OBResidueIter_SetName;
*SetNum = *Chemistry::OpenBabelc::OBResidueIter_SetNum;
*SetChain = *Chemistry::OpenBabelc::OBResidueIter_SetChain;
*SetChainNum = *Chemistry::OpenBabelc::OBResidueIter_SetChainNum;
*SetIdx = *Chemistry::OpenBabelc::OBResidueIter_SetIdx;
*SetAtomID = *Chemistry::OpenBabelc::OBResidueIter_SetAtomID;
*SetHetAtom = *Chemistry::OpenBabelc::OBResidueIter_SetHetAtom;
*SetSerialNum = *Chemistry::OpenBabelc::OBResidueIter_SetSerialNum;
*GetName = *Chemistry::OpenBabelc::OBResidueIter_GetName;
*GetNum = *Chemistry::OpenBabelc::OBResidueIter_GetNum;
*GetNumAtoms = *Chemistry::OpenBabelc::OBResidueIter_GetNumAtoms;
*GetChain = *Chemistry::OpenBabelc::OBResidueIter_GetChain;
*GetChainNum = *Chemistry::OpenBabelc::OBResidueIter_GetChainNum;
*GetIdx = *Chemistry::OpenBabelc::OBResidueIter_GetIdx;
*GetResKey = *Chemistry::OpenBabelc::OBResidueIter_GetResKey;
*GetAtoms = *Chemistry::OpenBabelc::OBResidueIter_GetAtoms;
*GetBonds = *Chemistry::OpenBabelc::OBResidueIter_GetBonds;
*GetAtomID = *Chemistry::OpenBabelc::OBResidueIter_GetAtomID;
*GetSerialNum = *Chemistry::OpenBabelc::OBResidueIter_GetSerialNum;
*GetAminoAcidProperty = *Chemistry::OpenBabelc::OBResidueIter_GetAminoAcidProperty;
*GetAtomProperty = *Chemistry::OpenBabelc::OBResidueIter_GetAtomProperty;
*GetResidueProperty = *Chemistry::OpenBabelc::OBResidueIter_GetResidueProperty;
*IsHetAtom = *Chemistry::OpenBabelc::OBResidueIter_IsHetAtom;
*IsResidueType = *Chemistry::OpenBabelc::OBResidueIter_IsResidueType;
*BeginAtoms = *Chemistry::OpenBabelc::OBResidueIter_BeginAtoms;
*EndAtoms = *Chemistry::OpenBabelc::OBResidueIter_EndAtoms;
*BeginAtom = *Chemistry::OpenBabelc::OBResidueIter_BeginAtom;
*NextAtom = *Chemistry::OpenBabelc::OBResidueIter_NextAtom;
*DoTransformations = *Chemistry::OpenBabelc::OBResidueIter_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBResidueIter_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBResidueIter_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBResidueIter_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBResidueIter_SetData;
*DataSize = *Chemistry::OpenBabelc::OBResidueIter_DataSize;
*GetData = *Chemistry::OpenBabelc::OBResidueIter_GetData;
*BeginData = *Chemistry::OpenBabelc::OBResidueIter_BeginData;
*EndData = *Chemistry::OpenBabelc::OBResidueIter_EndData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBResidueAtomIter ##############

package Chemistry::OpenBabel::OBResidueAtomIter;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBResidueAtomIter(@_);
    bless $self, $pkg if defined($self);
}

*good = *Chemistry::OpenBabelc::OBResidueAtomIter_good;
*inc = *Chemistry::OpenBabelc::OBResidueAtomIter_inc;
*deref = *Chemistry::OpenBabelc::OBResidueAtomIter_deref;
*__ref__ = *Chemistry::OpenBabelc::OBResidueAtomIter___ref__;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBResidueAtomIter($self);
        delete $OWNER{$self};
    }
}

*swig_Visit_get = *Chemistry::OpenBabelc::OBResidueAtomIter_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBResidueAtomIter_Visit_set;
*Clear = *Chemistry::OpenBabelc::OBResidueAtomIter_Clear;
*SetIdx = *Chemistry::OpenBabelc::OBResidueAtomIter_SetIdx;
*SetHyb = *Chemistry::OpenBabelc::OBResidueAtomIter_SetHyb;
*SetAtomicNum = *Chemistry::OpenBabelc::OBResidueAtomIter_SetAtomicNum;
*SetIsotope = *Chemistry::OpenBabelc::OBResidueAtomIter_SetIsotope;
*SetImplicitValence = *Chemistry::OpenBabelc::OBResidueAtomIter_SetImplicitValence;
*IncrementImplicitValence = *Chemistry::OpenBabelc::OBResidueAtomIter_IncrementImplicitValence;
*DecrementImplicitValence = *Chemistry::OpenBabelc::OBResidueAtomIter_DecrementImplicitValence;
*SetFormalCharge = *Chemistry::OpenBabelc::OBResidueAtomIter_SetFormalCharge;
*SetSpinMultiplicity = *Chemistry::OpenBabelc::OBResidueAtomIter_SetSpinMultiplicity;
*SetType = *Chemistry::OpenBabelc::OBResidueAtomIter_SetType;
*SetPartialCharge = *Chemistry::OpenBabelc::OBResidueAtomIter_SetPartialCharge;
*SetVector = *Chemistry::OpenBabelc::OBResidueAtomIter_SetVector;
*SetCoordPtr = *Chemistry::OpenBabelc::OBResidueAtomIter_SetCoordPtr;
*SetResidue = *Chemistry::OpenBabelc::OBResidueAtomIter_SetResidue;
*SetParent = *Chemistry::OpenBabelc::OBResidueAtomIter_SetParent;
*SetAromatic = *Chemistry::OpenBabelc::OBResidueAtomIter_SetAromatic;
*UnsetAromatic = *Chemistry::OpenBabelc::OBResidueAtomIter_UnsetAromatic;
*SetClockwiseStereo = *Chemistry::OpenBabelc::OBResidueAtomIter_SetClockwiseStereo;
*SetAntiClockwiseStereo = *Chemistry::OpenBabelc::OBResidueAtomIter_SetAntiClockwiseStereo;
*SetPositiveStereo = *Chemistry::OpenBabelc::OBResidueAtomIter_SetPositiveStereo;
*SetNegativeStereo = *Chemistry::OpenBabelc::OBResidueAtomIter_SetNegativeStereo;
*UnsetStereo = *Chemistry::OpenBabelc::OBResidueAtomIter_UnsetStereo;
*SetInRing = *Chemistry::OpenBabelc::OBResidueAtomIter_SetInRing;
*SetChiral = *Chemistry::OpenBabelc::OBResidueAtomIter_SetChiral;
*ClearCoordPtr = *Chemistry::OpenBabelc::OBResidueAtomIter_ClearCoordPtr;
*GetFormalCharge = *Chemistry::OpenBabelc::OBResidueAtomIter_GetFormalCharge;
*GetAtomicNum = *Chemistry::OpenBabelc::OBResidueAtomIter_GetAtomicNum;
*GetIsotope = *Chemistry::OpenBabelc::OBResidueAtomIter_GetIsotope;
*GetSpinMultiplicity = *Chemistry::OpenBabelc::OBResidueAtomIter_GetSpinMultiplicity;
*GetAtomicMass = *Chemistry::OpenBabelc::OBResidueAtomIter_GetAtomicMass;
*GetExactMass = *Chemistry::OpenBabelc::OBResidueAtomIter_GetExactMass;
*GetIdx = *Chemistry::OpenBabelc::OBResidueAtomIter_GetIdx;
*GetCoordinateIdx = *Chemistry::OpenBabelc::OBResidueAtomIter_GetCoordinateIdx;
*GetCIdx = *Chemistry::OpenBabelc::OBResidueAtomIter_GetCIdx;
*GetValence = *Chemistry::OpenBabelc::OBResidueAtomIter_GetValence;
*GetHyb = *Chemistry::OpenBabelc::OBResidueAtomIter_GetHyb;
*GetImplicitValence = *Chemistry::OpenBabelc::OBResidueAtomIter_GetImplicitValence;
*GetHvyValence = *Chemistry::OpenBabelc::OBResidueAtomIter_GetHvyValence;
*GetHeteroValence = *Chemistry::OpenBabelc::OBResidueAtomIter_GetHeteroValence;
*GetType = *Chemistry::OpenBabelc::OBResidueAtomIter_GetType;
*GetX = *Chemistry::OpenBabelc::OBResidueAtomIter_GetX;
*x = *Chemistry::OpenBabelc::OBResidueAtomIter_x;
*GetY = *Chemistry::OpenBabelc::OBResidueAtomIter_GetY;
*y = *Chemistry::OpenBabelc::OBResidueAtomIter_y;
*GetZ = *Chemistry::OpenBabelc::OBResidueAtomIter_GetZ;
*z = *Chemistry::OpenBabelc::OBResidueAtomIter_z;
*GetCoordinate = *Chemistry::OpenBabelc::OBResidueAtomIter_GetCoordinate;
*GetVector = *Chemistry::OpenBabelc::OBResidueAtomIter_GetVector;
*GetPartialCharge = *Chemistry::OpenBabelc::OBResidueAtomIter_GetPartialCharge;
*GetResidue = *Chemistry::OpenBabelc::OBResidueAtomIter_GetResidue;
*GetParent = *Chemistry::OpenBabelc::OBResidueAtomIter_GetParent;
*GetNewBondVector = *Chemistry::OpenBabelc::OBResidueAtomIter_GetNewBondVector;
*GetBond = *Chemistry::OpenBabelc::OBResidueAtomIter_GetBond;
*GetNextAtom = *Chemistry::OpenBabelc::OBResidueAtomIter_GetNextAtom;
*BeginBonds = *Chemistry::OpenBabelc::OBResidueAtomIter_BeginBonds;
*EndBonds = *Chemistry::OpenBabelc::OBResidueAtomIter_EndBonds;
*BeginBond = *Chemistry::OpenBabelc::OBResidueAtomIter_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBResidueAtomIter_NextBond;
*BeginNbrAtom = *Chemistry::OpenBabelc::OBResidueAtomIter_BeginNbrAtom;
*NextNbrAtom = *Chemistry::OpenBabelc::OBResidueAtomIter_NextNbrAtom;
*GetDistance = *Chemistry::OpenBabelc::OBResidueAtomIter_GetDistance;
*GetAngle = *Chemistry::OpenBabelc::OBResidueAtomIter_GetAngle;
*NewResidue = *Chemistry::OpenBabelc::OBResidueAtomIter_NewResidue;
*DeleteResidue = *Chemistry::OpenBabelc::OBResidueAtomIter_DeleteResidue;
*AddBond = *Chemistry::OpenBabelc::OBResidueAtomIter_AddBond;
*InsertBond = *Chemistry::OpenBabelc::OBResidueAtomIter_InsertBond;
*DeleteBond = *Chemistry::OpenBabelc::OBResidueAtomIter_DeleteBond;
*ClearBond = *Chemistry::OpenBabelc::OBResidueAtomIter_ClearBond;
*CountFreeOxygens = *Chemistry::OpenBabelc::OBResidueAtomIter_CountFreeOxygens;
*ImplicitHydrogenCount = *Chemistry::OpenBabelc::OBResidueAtomIter_ImplicitHydrogenCount;
*ExplicitHydrogenCount = *Chemistry::OpenBabelc::OBResidueAtomIter_ExplicitHydrogenCount;
*MemberOfRingCount = *Chemistry::OpenBabelc::OBResidueAtomIter_MemberOfRingCount;
*MemberOfRingSize = *Chemistry::OpenBabelc::OBResidueAtomIter_MemberOfRingSize;
*CountRingBonds = *Chemistry::OpenBabelc::OBResidueAtomIter_CountRingBonds;
*SmallestBondAngle = *Chemistry::OpenBabelc::OBResidueAtomIter_SmallestBondAngle;
*AverageBondAngle = *Chemistry::OpenBabelc::OBResidueAtomIter_AverageBondAngle;
*BOSum = *Chemistry::OpenBabelc::OBResidueAtomIter_BOSum;
*KBOSum = *Chemistry::OpenBabelc::OBResidueAtomIter_KBOSum;
*HtoMethyl = *Chemistry::OpenBabelc::OBResidueAtomIter_HtoMethyl;
*SetHybAndGeom = *Chemistry::OpenBabelc::OBResidueAtomIter_SetHybAndGeom;
*ForceNoH = *Chemistry::OpenBabelc::OBResidueAtomIter_ForceNoH;
*HasNoHForced = *Chemistry::OpenBabelc::OBResidueAtomIter_HasNoHForced;
*HasResidue = *Chemistry::OpenBabelc::OBResidueAtomIter_HasResidue;
*IsHydrogen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsHydrogen;
*IsCarbon = *Chemistry::OpenBabelc::OBResidueAtomIter_IsCarbon;
*IsNitrogen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsNitrogen;
*IsOxygen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsOxygen;
*IsSulfur = *Chemistry::OpenBabelc::OBResidueAtomIter_IsSulfur;
*IsPhosphorus = *Chemistry::OpenBabelc::OBResidueAtomIter_IsPhosphorus;
*IsAromatic = *Chemistry::OpenBabelc::OBResidueAtomIter_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBResidueAtomIter_IsInRing;
*IsInRingSize = *Chemistry::OpenBabelc::OBResidueAtomIter_IsInRingSize;
*IsHeteroatom = *Chemistry::OpenBabelc::OBResidueAtomIter_IsHeteroatom;
*IsNotCorH = *Chemistry::OpenBabelc::OBResidueAtomIter_IsNotCorH;
*IsConnected = *Chemistry::OpenBabelc::OBResidueAtomIter_IsConnected;
*IsOneThree = *Chemistry::OpenBabelc::OBResidueAtomIter_IsOneThree;
*IsOneFour = *Chemistry::OpenBabelc::OBResidueAtomIter_IsOneFour;
*IsCarboxylOxygen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsCarboxylOxygen;
*IsPhosphateOxygen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsPhosphateOxygen;
*IsSulfateOxygen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsSulfateOxygen;
*IsNitroOxygen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsNitroOxygen;
*IsAmideNitrogen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsAmideNitrogen;
*IsPolarHydrogen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsPolarHydrogen;
*IsNonPolarHydrogen = *Chemistry::OpenBabelc::OBResidueAtomIter_IsNonPolarHydrogen;
*IsAromaticNOxide = *Chemistry::OpenBabelc::OBResidueAtomIter_IsAromaticNOxide;
*IsChiral = *Chemistry::OpenBabelc::OBResidueAtomIter_IsChiral;
*IsAxial = *Chemistry::OpenBabelc::OBResidueAtomIter_IsAxial;
*IsClockwise = *Chemistry::OpenBabelc::OBResidueAtomIter_IsClockwise;
*IsAntiClockwise = *Chemistry::OpenBabelc::OBResidueAtomIter_IsAntiClockwise;
*IsPositiveStereo = *Chemistry::OpenBabelc::OBResidueAtomIter_IsPositiveStereo;
*IsNegativeStereo = *Chemistry::OpenBabelc::OBResidueAtomIter_IsNegativeStereo;
*HasChiralitySpecified = *Chemistry::OpenBabelc::OBResidueAtomIter_HasChiralitySpecified;
*HasChiralVolume = *Chemistry::OpenBabelc::OBResidueAtomIter_HasChiralVolume;
*IsHbondAcceptor = *Chemistry::OpenBabelc::OBResidueAtomIter_IsHbondAcceptor;
*IsHbondDonor = *Chemistry::OpenBabelc::OBResidueAtomIter_IsHbondDonor;
*IsHbondDonorH = *Chemistry::OpenBabelc::OBResidueAtomIter_IsHbondDonorH;
*HasAlphaBetaUnsat = *Chemistry::OpenBabelc::OBResidueAtomIter_HasAlphaBetaUnsat;
*HasBondOfOrder = *Chemistry::OpenBabelc::OBResidueAtomIter_HasBondOfOrder;
*CountBondsOfOrder = *Chemistry::OpenBabelc::OBResidueAtomIter_CountBondsOfOrder;
*HasNonSingleBond = *Chemistry::OpenBabelc::OBResidueAtomIter_HasNonSingleBond;
*HasSingleBond = *Chemistry::OpenBabelc::OBResidueAtomIter_HasSingleBond;
*HasDoubleBond = *Chemistry::OpenBabelc::OBResidueAtomIter_HasDoubleBond;
*HasAromaticBond = *Chemistry::OpenBabelc::OBResidueAtomIter_HasAromaticBond;
*MatchesSMARTS = *Chemistry::OpenBabelc::OBResidueAtomIter_MatchesSMARTS;
*DoTransformations = *Chemistry::OpenBabelc::OBResidueAtomIter_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBResidueAtomIter_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBResidueAtomIter_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBResidueAtomIter_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBResidueAtomIter_SetData;
*DataSize = *Chemistry::OpenBabelc::OBResidueAtomIter_DataSize;
*GetData = *Chemistry::OpenBabelc::OBResidueAtomIter_GetData;
*BeginData = *Chemistry::OpenBabelc::OBResidueAtomIter_BeginData;
*EndData = *Chemistry::OpenBabelc::OBResidueAtomIter_EndData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


# ------- VARIABLE STUBS --------

package Chemistry::OpenBabel;

*FILE_SEP_CHAR = *Chemistry::OpenBabelc::FILE_SEP_CHAR;

my %__VZero_hash;
tie %__VZero_hash,"Chemistry::OpenBabel::vector3", $Chemistry::OpenBabelc::VZero;
$VZero= \%__VZero_hash;
bless $VZero, Chemistry::OpenBabel::vector3;

my %__VX_hash;
tie %__VX_hash,"Chemistry::OpenBabel::vector3", $Chemistry::OpenBabelc::VX;
$VX= \%__VX_hash;
bless $VX, Chemistry::OpenBabel::vector3;

my %__VY_hash;
tie %__VY_hash,"Chemistry::OpenBabel::vector3", $Chemistry::OpenBabelc::VY;
$VY= \%__VY_hash;
bless $VY, Chemistry::OpenBabel::vector3;

my %__VZ_hash;
tie %__VZ_hash,"Chemistry::OpenBabel::vector3", $Chemistry::OpenBabelc::VZ;
$VZ= \%__VZ_hash;
bless $VZ, Chemistry::OpenBabel::vector3;
*UndefinedData = *Chemistry::OpenBabelc::UndefinedData;
*PairData = *Chemistry::OpenBabelc::PairData;
*EnergyData = *Chemistry::OpenBabelc::EnergyData;
*CommentData = *Chemistry::OpenBabelc::CommentData;
*ConformerData = *Chemistry::OpenBabelc::ConformerData;
*ExternalBondData = *Chemistry::OpenBabelc::ExternalBondData;
*RotamerList = *Chemistry::OpenBabelc::RotamerList;
*VirtualBondData = *Chemistry::OpenBabelc::VirtualBondData;
*RingData = *Chemistry::OpenBabelc::RingData;
*TorsionData = *Chemistry::OpenBabelc::TorsionData;
*AngleData = *Chemistry::OpenBabelc::AngleData;
*SerialNums = *Chemistry::OpenBabelc::SerialNums;
*UnitCell = *Chemistry::OpenBabelc::UnitCell;
*SpinData = *Chemistry::OpenBabelc::SpinData;
*ChargeData = *Chemistry::OpenBabelc::ChargeData;
*SymmetryData = *Chemistry::OpenBabelc::SymmetryData;
*ChiralData = *Chemistry::OpenBabelc::ChiralData;
*OccupationData = *Chemistry::OpenBabelc::OccupationData;
*DensityData = *Chemistry::OpenBabelc::DensityData;
*ElectronicData = *Chemistry::OpenBabelc::ElectronicData;
*VibrationData = *Chemistry::OpenBabelc::VibrationData;
*RotationData = *Chemistry::OpenBabelc::RotationData;
*NuclearData = *Chemistry::OpenBabelc::NuclearData;
*SetData = *Chemistry::OpenBabelc::SetData;
*CustomData0 = *Chemistry::OpenBabelc::CustomData0;
*CustomData1 = *Chemistry::OpenBabelc::CustomData1;
*CustomData2 = *Chemistry::OpenBabelc::CustomData2;
*CustomData3 = *Chemistry::OpenBabelc::CustomData3;
*CustomData4 = *Chemistry::OpenBabelc::CustomData4;
*CustomData5 = *Chemistry::OpenBabelc::CustomData5;
*CustomData6 = *Chemistry::OpenBabelc::CustomData6;
*CustomData7 = *Chemistry::OpenBabelc::CustomData7;
*CustomData8 = *Chemistry::OpenBabelc::CustomData8;
*CustomData9 = *Chemistry::OpenBabelc::CustomData9;
*CustomData10 = *Chemistry::OpenBabelc::CustomData10;
*CustomData11 = *Chemistry::OpenBabelc::CustomData11;
*CustomData12 = *Chemistry::OpenBabelc::CustomData12;
*CustomData13 = *Chemistry::OpenBabelc::CustomData13;
*CustomData14 = *Chemistry::OpenBabelc::CustomData14;
*CustomData15 = *Chemistry::OpenBabelc::CustomData15;
*output = *Chemistry::OpenBabelc::output;
*input = *Chemistry::OpenBabelc::input;
*calcvolume = *Chemistry::OpenBabelc::calcvolume;
*obError = *Chemistry::OpenBabelc::obError;
*obWarning = *Chemistry::OpenBabelc::obWarning;
*obInfo = *Chemistry::OpenBabelc::obInfo;
*obAuditMsg = *Chemistry::OpenBabelc::obAuditMsg;
*obDebug = *Chemistry::OpenBabelc::obDebug;

my %__obErrorLog_hash;
tie %__obErrorLog_hash,"Chemistry::OpenBabel::OBMessageHandler", $Chemistry::OpenBabelc::obErrorLog;
$obErrorLog= \%__obErrorLog_hash;
bless $obErrorLog, Chemistry::OpenBabel::OBMessageHandler;
*NOTREADABLE = *Chemistry::OpenBabelc::NOTREADABLE;
*READONEONLY = *Chemistry::OpenBabelc::READONEONLY;
*READBINARY = *Chemistry::OpenBabelc::READBINARY;
*ZEROATOMSOK = *Chemistry::OpenBabelc::ZEROATOMSOK;
*NOTWRITABLE = *Chemistry::OpenBabelc::NOTWRITABLE;
*WRITEONEONLY = *Chemistry::OpenBabelc::WRITEONEONLY;
*WRITEBINARY = *Chemistry::OpenBabelc::WRITEBINARY;
*DEFAULTFORMAT = *Chemistry::OpenBabelc::DEFAULTFORMAT;
*MAXSETNO = *Chemistry::OpenBabelc::MAXSETNO;
*MAXELEM = *Chemistry::OpenBabelc::MAXELEM;
*MINELEM = *Chemistry::OpenBabelc::MINELEM;
*MAXRES = *Chemistry::OpenBabelc::MAXRES;
*MINRES = *Chemistry::OpenBabelc::MINRES;
*AA_ALA = *Chemistry::OpenBabelc::AA_ALA;
*AA_GLY = *Chemistry::OpenBabelc::AA_GLY;
*AA_LEU = *Chemistry::OpenBabelc::AA_LEU;
*AA_SER = *Chemistry::OpenBabelc::AA_SER;
*AA_VAL = *Chemistry::OpenBabelc::AA_VAL;
*AA_THR = *Chemistry::OpenBabelc::AA_THR;
*AA_LYS = *Chemistry::OpenBabelc::AA_LYS;
*AA_ASP = *Chemistry::OpenBabelc::AA_ASP;
*AA_ILE = *Chemistry::OpenBabelc::AA_ILE;
*AA_ASN = *Chemistry::OpenBabelc::AA_ASN;
*AA_GLU = *Chemistry::OpenBabelc::AA_GLU;
*AA_PRO = *Chemistry::OpenBabelc::AA_PRO;
*AA_ARG = *Chemistry::OpenBabelc::AA_ARG;
*AA_PHE = *Chemistry::OpenBabelc::AA_PHE;
*AA_GLN = *Chemistry::OpenBabelc::AA_GLN;
*AA_TYR = *Chemistry::OpenBabelc::AA_TYR;
*AA_HIS = *Chemistry::OpenBabelc::AA_HIS;
*AA_CYS = *Chemistry::OpenBabelc::AA_CYS;
*AA_MET = *Chemistry::OpenBabelc::AA_MET;
*AA_TRP = *Chemistry::OpenBabelc::AA_TRP;
*ACIDIC = *Chemistry::OpenBabelc::ACIDIC;
*ACYCLIC = *Chemistry::OpenBabelc::ACYCLIC;
*ALIPHATIC = *Chemistry::OpenBabelc::ALIPHATIC;
*AROMATIC = *Chemistry::OpenBabelc::AROMATIC;
*BASIC = *Chemistry::OpenBabelc::BASIC;
*BURIED = *Chemistry::OpenBabelc::BURIED;
*CHARGED = *Chemistry::OpenBabelc::CHARGED;
*CYCLIC = *Chemistry::OpenBabelc::CYCLIC;
*HYDROPHOBIC = *Chemistry::OpenBabelc::HYDROPHOBIC;
*LARGE = *Chemistry::OpenBabelc::LARGE;
*MEDIUM = *Chemistry::OpenBabelc::MEDIUM;
*NEGATIVE = *Chemistry::OpenBabelc::NEGATIVE;
*NEUTRAL = *Chemistry::OpenBabelc::NEUTRAL;
*POLAR = *Chemistry::OpenBabelc::POLAR;
*POSITIVE = *Chemistry::OpenBabelc::POSITIVE;
*SMALL = *Chemistry::OpenBabelc::SMALL;
*SURFACE = *Chemistry::OpenBabelc::SURFACE;
*ALPHA_CARBON = *Chemistry::OpenBabelc::ALPHA_CARBON;
*AMINO_BACKBONE = *Chemistry::OpenBabelc::AMINO_BACKBONE;
*BACKBONE = *Chemistry::OpenBabelc::BACKBONE;
*CYSTEINE_SULPHUR = *Chemistry::OpenBabelc::CYSTEINE_SULPHUR;
*LIGAND = *Chemistry::OpenBabelc::LIGAND;
*NUCLEIC_BACKBONE = *Chemistry::OpenBabelc::NUCLEIC_BACKBONE;
*SHAPELY_BACKBONE = *Chemistry::OpenBabelc::SHAPELY_BACKBONE;
*SHAPELY_SPECIAL = *Chemistry::OpenBabelc::SHAPELY_SPECIAL;
*SIDECHAIN = *Chemistry::OpenBabelc::SIDECHAIN;
*SUGAR_PHOSPHATE = *Chemistry::OpenBabelc::SUGAR_PHOSPHATE;
*ALA = *Chemistry::OpenBabelc::ALA;
*GLY = *Chemistry::OpenBabelc::GLY;
*LEU = *Chemistry::OpenBabelc::LEU;
*SER = *Chemistry::OpenBabelc::SER;
*VAL = *Chemistry::OpenBabelc::VAL;
*THR = *Chemistry::OpenBabelc::THR;
*LYS = *Chemistry::OpenBabelc::LYS;
*ASP = *Chemistry::OpenBabelc::ASP;
*ILE = *Chemistry::OpenBabelc::ILE;
*ASN = *Chemistry::OpenBabelc::ASN;
*GLU = *Chemistry::OpenBabelc::GLU;
*PRO = *Chemistry::OpenBabelc::PRO;
*ARG = *Chemistry::OpenBabelc::ARG;
*PHE = *Chemistry::OpenBabelc::PHE;
*GLN = *Chemistry::OpenBabelc::GLN;
*TYR = *Chemistry::OpenBabelc::TYR;
*HIS = *Chemistry::OpenBabelc::HIS;
*CYS = *Chemistry::OpenBabelc::CYS;
*MET = *Chemistry::OpenBabelc::MET;
*TRP = *Chemistry::OpenBabelc::TRP;
*ASX = *Chemistry::OpenBabelc::ASX;
*GLX = *Chemistry::OpenBabelc::GLX;
*PCA = *Chemistry::OpenBabelc::PCA;
*HYP = *Chemistry::OpenBabelc::HYP;
*A = *Chemistry::OpenBabelc::A;
*C = *Chemistry::OpenBabelc::C;
*G = *Chemistry::OpenBabelc::G;
*T = *Chemistry::OpenBabelc::T;
*U = *Chemistry::OpenBabelc::U;
*UPLUS = *Chemistry::OpenBabelc::UPLUS;
*I = *Chemistry::OpenBabelc::I;
*OMC = *Chemistry::OpenBabelc::OMC;
*M2G = *Chemistry::OpenBabelc::M2G;
*OMG = *Chemistry::OpenBabelc::OMG;
*YG = *Chemistry::OpenBabelc::YG;
*H2U = *Chemistry::OpenBabelc::H2U;
*PSU = *Chemistry::OpenBabelc::PSU;
*UNK = *Chemistry::OpenBabelc::UNK;
*ACE = *Chemistry::OpenBabelc::ACE;
*FOR = *Chemistry::OpenBabelc::FOR;
*HOH = *Chemistry::OpenBabelc::HOH;
*DOD = *Chemistry::OpenBabelc::DOD;
*SO4 = *Chemistry::OpenBabelc::SO4;
*PO4 = *Chemistry::OpenBabelc::PO4;
*NAD = *Chemistry::OpenBabelc::NAD;
*COA = *Chemistry::OpenBabelc::COA;
*NAP = *Chemistry::OpenBabelc::NAP;
*NDP = *Chemistry::OpenBabelc::NDP;
*AMINO = *Chemistry::OpenBabelc::AMINO;
*AMINO_NUCLEO = *Chemistry::OpenBabelc::AMINO_NUCLEO;
*COENZYME = *Chemistry::OpenBabelc::COENZYME;
*ION = *Chemistry::OpenBabelc::ION;
*NUCLEO = *Chemistry::OpenBabelc::NUCLEO;
*PROTEIN = *Chemistry::OpenBabelc::PROTEIN;
*PURINE = *Chemistry::OpenBabelc::PURINE;
*PYRIMIDINE = *Chemistry::OpenBabelc::PYRIMIDINE;
*SOLVENT = *Chemistry::OpenBabelc::SOLVENT;
*WATER = *Chemistry::OpenBabelc::WATER;
*Residue = *Chemistry::OpenBabelc::Residue;
*ElemDesc = *Chemistry::OpenBabelc::ElemDesc;
*ResNo = *Chemistry::OpenBabelc::ResNo;
*ElemNo = *Chemistry::OpenBabelc::ElemNo;
*OB_4RING_ATOM = *Chemistry::OpenBabelc::OB_4RING_ATOM;
*OB_3RING_ATOM = *Chemistry::OpenBabelc::OB_3RING_ATOM;
*OB_AROMATIC_ATOM = *Chemistry::OpenBabelc::OB_AROMATIC_ATOM;
*OB_RING_ATOM = *Chemistry::OpenBabelc::OB_RING_ATOM;
*OB_CSTEREO_ATOM = *Chemistry::OpenBabelc::OB_CSTEREO_ATOM;
*OB_ACSTEREO_ATOM = *Chemistry::OpenBabelc::OB_ACSTEREO_ATOM;
*OB_DONOR_ATOM = *Chemistry::OpenBabelc::OB_DONOR_ATOM;
*OB_ACCEPTOR_ATOM = *Chemistry::OpenBabelc::OB_ACCEPTOR_ATOM;
*OB_CHIRAL_ATOM = *Chemistry::OpenBabelc::OB_CHIRAL_ATOM;
*OB_POS_CHIRAL_ATOM = *Chemistry::OpenBabelc::OB_POS_CHIRAL_ATOM;
*OB_NEG_CHIRAL_ATOM = *Chemistry::OpenBabelc::OB_NEG_CHIRAL_ATOM;
*OB_ATOM_HAS_NO_H = *Chemistry::OpenBabelc::OB_ATOM_HAS_NO_H;
*OB_AROMATIC_BOND = *Chemistry::OpenBabelc::OB_AROMATIC_BOND;
*OB_WEDGE_BOND = *Chemistry::OpenBabelc::OB_WEDGE_BOND;
*OB_HASH_BOND = *Chemistry::OpenBabelc::OB_HASH_BOND;
*OB_RING_BOND = *Chemistry::OpenBabelc::OB_RING_BOND;
*OB_TORUP_BOND = *Chemistry::OpenBabelc::OB_TORUP_BOND;
*OB_TORDOWN_BOND = *Chemistry::OpenBabelc::OB_TORDOWN_BOND;
*OB_KSINGLE_BOND = *Chemistry::OpenBabelc::OB_KSINGLE_BOND;
*OB_KDOUBLE_BOND = *Chemistry::OpenBabelc::OB_KDOUBLE_BOND;
*OB_KTRIPLE_BOND = *Chemistry::OpenBabelc::OB_KTRIPLE_BOND;
*OB_CLOSURE_BOND = *Chemistry::OpenBabelc::OB_CLOSURE_BOND;
*OB_SSSR_MOL = *Chemistry::OpenBabelc::OB_SSSR_MOL;
*OB_RINGFLAGS_MOL = *Chemistry::OpenBabelc::OB_RINGFLAGS_MOL;
*OB_AROMATIC_MOL = *Chemistry::OpenBabelc::OB_AROMATIC_MOL;
*OB_ATOMTYPES_MOL = *Chemistry::OpenBabelc::OB_ATOMTYPES_MOL;
*OB_CHIRALITY_MOL = *Chemistry::OpenBabelc::OB_CHIRALITY_MOL;
*OB_PCHARGE_MOL = *Chemistry::OpenBabelc::OB_PCHARGE_MOL;
*OB_HYBRID_MOL = *Chemistry::OpenBabelc::OB_HYBRID_MOL;
*OB_IMPVAL_MOL = *Chemistry::OpenBabelc::OB_IMPVAL_MOL;
*OB_KEKULE_MOL = *Chemistry::OpenBabelc::OB_KEKULE_MOL;
*OB_CLOSURE_MOL = *Chemistry::OpenBabelc::OB_CLOSURE_MOL;
*OB_H_ADDED_MOL = *Chemistry::OpenBabelc::OB_H_ADDED_MOL;
*OB_PH_CORRECTED_MOL = *Chemistry::OpenBabelc::OB_PH_CORRECTED_MOL;
*OB_AROM_CORRECTED_MOL = *Chemistry::OpenBabelc::OB_AROM_CORRECTED_MOL;
*OB_CHAINS_MOL = *Chemistry::OpenBabelc::OB_CHAINS_MOL;
*OB_TCHARGE_MOL = *Chemistry::OpenBabelc::OB_TCHARGE_MOL;
*OB_TSPIN_MOL = *Chemistry::OpenBabelc::OB_TSPIN_MOL;
*OB_CURRENT_CONFORMER = *Chemistry::OpenBabelc::OB_CURRENT_CONFORMER;

my %__etab_hash;
tie %__etab_hash,"Chemistry::OpenBabel::OBElementTable", $Chemistry::OpenBabelc::etab;
$etab= \%__etab_hash;
bless $etab, Chemistry::OpenBabel::OBElementTable;

my %__ttab_hash;
tie %__ttab_hash,"Chemistry::OpenBabel::OBTypeTable", $Chemistry::OpenBabelc::ttab;
$ttab= \%__ttab_hash;
bless $ttab, Chemistry::OpenBabel::OBTypeTable;

my %__isotab_hash;
tie %__isotab_hash,"Chemistry::OpenBabel::OBIsotopeTable", $Chemistry::OpenBabelc::isotab;
$isotab= \%__isotab_hash;
bless $isotab, Chemistry::OpenBabel::OBIsotopeTable;
*aromtyper = *Chemistry::OpenBabelc::aromtyper;
*atomtyper = *Chemistry::OpenBabelc::atomtyper;
*chainsparser = *Chemistry::OpenBabelc::chainsparser;

my %__resdat_hash;
tie %__resdat_hash,"Chemistry::OpenBabel::OBResidueData", $Chemistry::OpenBabelc::resdat;
$resdat= \%__resdat_hash;
bless $resdat, Chemistry::OpenBabel::OBResidueData;
*BUFF_SIZE = *Chemistry::OpenBabelc::BUFF_SIZE;
1;
