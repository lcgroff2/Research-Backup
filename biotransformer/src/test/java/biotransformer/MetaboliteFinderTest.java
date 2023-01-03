package biotransformer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import biotransformer.fingerprint.ChemStructureFingerprinter;
import biotransformer.utils.MetaboliteFinder;
import biotransformer.utils.MetaboliteFinder.FinderOption;

public class MetaboliteFinderTest extends MetaboliteFinder {

	
	public MetaboliteFinderTest() throws IOException, ParseException, CDKException {
		// TODO Auto-generated constructor stub
		super();
	}
	
	public static void main(String[] args) throws Exception{
		MetaboliteFinderTest mft = new MetaboliteFinderTest();
		IAtomContainer atc = mft.hsbt.smiParser.parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1");
		ArrayList<String> masses = new ArrayList<String>();
		ArrayList<String> formulas = new ArrayList<String>();

		masses.add("304.0941");
		masses.add("342.0946");
		
		mft.findSuperbioMetabolites(atc, masses, 0.01, true, "data/epicatechin_metabolites_identification_by_mass.sdf", FinderOption.MASS);
	
		formulas.add("C11H12O3");
		mft.findAllHumanMetabolites(atc, formulas, 0.01, 2, false, "data/epicatechin_metabolites_identification_by_formula.sdf", FinderOption.FORMULA);
		
	}

}
