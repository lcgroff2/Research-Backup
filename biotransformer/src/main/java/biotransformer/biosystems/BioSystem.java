/**
 * This class implements the class of biosystems. They can represent either individual
 * species/organisms or a collection thereof
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */
package biotransformer.biosystems;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.JsonParser.Feature;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;
import org.openscience.cdk.silent.SilentChemObjectBuilder;


import ambit2.smarts.SMIRKSManager;
import biotransformer.biomolecule.Enzyme;
import biotransformer.biomolecule.Enzyme.EnzymeName;
import biotransformer.transformation.MRPatterns;
import biotransformer.transformation.MReactionsFilter;
import biotransformer.transformation.MetabolicPathway;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.transformation.MRPatterns.ReactionName;
import biotransformer.transformation.MetabolicPathway.MPathwayName;

/**
 * @author Yannick Djoumbou Feunang
 *
 */
public class BioSystem {

	/**
	 * 
	 */
//	protected LinkedHashMap<String, ArrayList<MetabolicReaction>>	reactions = new LinkedHashMap<String, ArrayList<MetabolicReaction>>();
	protected LinkedHashMap<String, MetabolicReaction> reactionsHash = new LinkedHashMap<String, MetabolicReaction> ();
	protected LinkedHashMap<String, Double>	reactionsOcurrenceRatios = new LinkedHashMap<String, Double>();
	protected ArrayList<Enzyme>	enzymes = new  ArrayList<Enzyme>();
	protected LinkedHashMap<EnzymeName, Enzyme>	enzymesHash = new  LinkedHashMap<EnzymeName, Enzyme>();
	
	// Update this to include references too
	protected LinkedHashMap<MPathwayName, ArrayList<Enzyme>> metPathwaysHash = new LinkedHashMap<MPathwayName,ArrayList<Enzyme>>();
	public MReactionsFilter mrFilter;
	protected SMIRKSManager 		smrkMan	= new SMIRKSManager(SilentChemObjectBuilder.getInstance());	
	public BioSystemName name;

	
	
//	public static void main (String[] args) throws JsonParseException, JsonMappingException, IOException{
//		ObjectMapper mapper = new ObjectMapper();
//		mapper.configure(Feature.ALLOW_COMMENTS, true);
//		mapper.configure(Feature.ALLOW_BACKSLASH_ESCAPING_ANY_CHARACTER, true);
////		mapper.configure(Feature.AUTO_CLOSE_SOURCE, true);
//		 Biosystem bios = new Biosystem(BioSystemName.HUMAN, mapper);	 
//		 System.out.println(bios.enzymes.get(0).getName());
//		 System.out.println(bios.reactionsHash.keySet().size());
//		 System.out.println(bios.metPathwaysHash);
//	}
	
	public BioSystem(BioSystemName bsName, ObjectMapper mapper) throws JsonParseException, JsonMappingException, IOException {
		this.name 			= bsName;
		this.mrFilter		= new MReactionsFilter(bsName);
		setEnzymesList(mapper);
//		setReactionsList();
	}
	
	
	public BioSystem(BioSystemName bsName, ArrayList<Enzyme> enzymes){
		this.name 			= bsName;
		this.enzymes		= enzymes;
//		this.mrFilter		= new MReactionsFilter(bsName);
//		setEnzymesList();
//		setReactionsList();
	}
	
	
	public enum BioSystemName {
		HUMAN, ENVMICRO, GUTMICRO
	}

	@SuppressWarnings("unchecked")
	private void setEnzymesList(ObjectMapper mapper) throws JsonParseException, JsonMappingException, IOException {
		
		/*
		 * Parsing the JSON files that constitute the database.
		 */				
//		System.err.println(BioSystem.class.getClassLoader().getSystemResource("enzymes.json"));
//		InputStream enzymes = BioSystem.class.getResourceAsStream("database/enzymes.json");
//		
//		ClassLoader classLoader = this.getClass().getClassLoader();
//		File enzymes = new File(Thread.currentThread().getContextClassLoader().getResource("database/enzymes.json").getFile());
//		File enzymes = new File("enzymes.json");
//		File biosystemEnzymes = new File("biosystemEnzymes.json");
//		File enzymeReactions = new File("enzymeReactions.json");
//		File metabolicReactions = new File("metabolicReactions.json");
//		File bioSystemsReactionORatios = new File("bioSystemsReactionORatios.json");
//		File pathways = new File("pathways.json");
//		File enzymes = new File(Thread.currentThread().getContextClassLoader().getResourceAsStream("file.txt").get);
		File enzymes = new File("database/enzymes.json");
		File biosystemEnzymes = new File("database/biosystemEnzymes.json");
		File enzymeReactions = new File("database/enzymeReactions.json");
		File metabolicReactions = new File("database/metabolicReactions.json");
		File bioSystemsReactionORatios = new File("database/biosystemsReactionORatios.json");
		File pathways = new File("database/pathways.json");
		
		// Parsing the list of enzymes and their attributes
		LinkedHashMap<String,Object> allEnzymes = (LinkedHashMap<String,Object>) mapper.readValue(enzymes, Map.class).get("enzymes");
//		System.out.println("Enzymes: " + allEnzymes);
		//		System.out.println("CYP1A2: " + allEnzymes.get("CYP1A2").get("description").toString());
		//		LinkedHashMap<String, Object> cyp1A2 = (LinkedHashMap<String, Object>) allEnzymes.get("CYP1A2");
		//		System.out.println("CYP1A2: " + cyp1A2.get("biosystems"));
		//		System.out.println("CYP1A2: " + allEnzymes.get("CYP1A2").get("biosystems").getClass());
		

		LinkedHashMap<String,Object> allRe = (LinkedHashMap<String, Object>) mapper.readValue(metabolicReactions, Map.class).get("reactions");
		//		System.out.println("REACTIONS: " + allRe);

		Map<String,Object> allEnzymesByBiosystem = mapper.readValue(biosystemEnzymes, Map.class);
		//		System.out.println("This organism: " + this.name.toString() );
		ArrayList<String> bioSysEnzymeList = (ArrayList<String>) ((LinkedHashMap<String, ArrayList<String>>) 
				allEnzymesByBiosystem.get("enzymeLists")).get(this.name.toString());
		//		System.out.println("Enzyme list for " + this.name.toString() + ": " + bioSysEnzymeList.size());
		
		Map<String,Object> allEnzToReactions = mapper.readValue(enzymeReactions, Map.class);
		LinkedHashMap<String, ArrayList<String>> enzymeReactionList = (LinkedHashMap<String, ArrayList<String>>) 
				allEnzToReactions.get("eReactionLists");

		//		System.out.println("enzymeReactionLists: " + enzymeReactionList);
		
		this.reactionsOcurrenceRatios = (LinkedHashMap<String, Double>) ((LinkedHashMap<String, Object>) 
				mapper.readValue(bioSystemsReactionORatios, Map.class).get("reactionsORatios")).get(this.name.toString());
		
		/*
		 * create a unique list of reactions and create them.
		 * This helps because a reaction can be catalyzed by many 
		 * enzymes, and we do not one to create the same object many times.
		 */		
		
			//		LinkedHashMap<String, MetabolicReaction> mReactionObjects = new LinkedHashMap<String, MetabolicReaction>();
		
		/*
		 * Now for each enzyme associated with the biosystem, built the reactions arraylist and then create the enzyme.
		 */		
		for( String e : bioSysEnzymeList ){			
			ArrayList<MetabolicReaction> enzymeSpecificReactionObjects = new ArrayList<MetabolicReaction>();
//						System.out.println(e);
			for(String s : enzymeReactionList.get(e)){
//								System.out.println(s);
				if(this.reactionsHash.containsKey(s)){
					enzymeSpecificReactionObjects.add(this.reactionsHash.get(s));
				} 				
				else {						
//					System.out.println(s);
					LinkedHashMap<String,Object> mrObj = (LinkedHashMap<String,Object>) allRe.get(s);
					String commonName = (String)mrObj.get("commonName");
					if(commonName == null || commonName.contentEquals("")){
						commonName = s;
					}
					String reactionBTMRID = (String)mrObj.get("btmrID");
//					System.err.println(reactionBTMRID);
					MetabolicReaction r = new MetabolicReaction(ReactionName.valueOf(s), commonName, reactionBTMRID, (String)mrObj.get("smirks"), 
							(ArrayList<String>)mrObj.get("smarts"), (ArrayList<String>)mrObj.get("negativeSmarts"), this.smrkMan);
					
					this.reactionsHash.put(s,r);
					enzymeSpecificReactionObjects.add(r);
//					r.display();
				}
			}
			
			 LinkedHashMap<String, Object> enz = (LinkedHashMap<String, Object>) allEnzymes.get(e);
//			 			 System.out.println(e);
			 LinkedHashMap<String, Object> biosystems = (LinkedHashMap<String, Object>) enz.get("biosystems");
//			 			 System.out.println(biosystems.get(this.name.toString()));
			 String description = (String) enz.get("description");
			 
			 
			 ArrayList<String> uniprot_ids = null;
			 ArrayList<String> cellularLocations = null;
			 if(biosystems.containsKey(this.name.toString())){
				 uniprot_ids = (ArrayList<String>)  ((LinkedHashMap<String, Object>) 
						 biosystems.get(this.name.toString())).get("uniprot_ids");
				 cellularLocations = (ArrayList<String>)  ((LinkedHashMap<String, Object>) 
						 biosystems.get(this.name.toString())).get("cellular_locations");
				 
				 //				 System.out.println(uniprot_ids);
				 //				 System.out.println(cellularLocations);				 
			 }

			 String acceptedName =  (String) enz.get("acceptedName");
			 if(acceptedName == null || acceptedName.contentEquals("")){
				 acceptedName =  (String)e;
			 }
			 
			Enzyme enzy = new Enzyme(e, description, uniprot_ids, cellularLocations, enzymeSpecificReactionObjects, acceptedName); 
			this.enzymes.add(enzy);
			this.enzymesHash.put(EnzymeName.valueOf(e), enzy);
//			enzy.display();
		}
		
//		System.out.println(reactionsHash);
//		System.err.println(reactionsHash.size());
		
//		System.err.println(this.reactionsHash);
		/*
		 * Now generating a pathway map.
		 */
		
		LinkedHashMap<String,Object> allPathways = (LinkedHashMap<String, Object>) mapper.readValue(pathways, Map.class).get("metabolicPathways");
		
		for(String pName : allPathways.keySet()) {
//						System.out.println(pName);
			LinkedHashMap<String,Object> lm = (LinkedHashMap<String,Object>) allPathways.get(pName);
//			System.out.println(lm.keySet().contains(this.name.toString()));
			
			if(lm.containsKey(this.name.toString())){
				LinkedHashMap<String,Object> bPaths = (LinkedHashMap<String,Object>) lm.get(this.name.toString());
				//				System.out.println(this.name.toString());
				//				System.out.println(lm);
				ArrayList<String> enz = (ArrayList<String>) bPaths.get("enzymes");
//				System.out.println(enz);
				
				this.metPathwaysHash.put(MPathwayName.valueOf(pName), new ArrayList<Enzyme>());
				
				for(String e_n : enz) {
					this.metPathwaysHash.get(MPathwayName.valueOf(pName)).add(enzymesHash.get(EnzymeName.valueOf(e_n)));
				}

			}

		}
		
//		System.out.println(metPathwaysHash);
		
		
	}

	public ArrayList<Enzyme> getEnzymesList(){
		return this.enzymes;
	}

	
	public SMIRKSManager getSmirksManager(){
		return this.smrkMan;
	}
	
	public  LinkedHashMap<String, MetabolicReaction> getReactionsHash(){
		return this.reactionsHash;
	}
	
	public LinkedHashMap<String, Double> getReactionsORatios(){
		return this.reactionsOcurrenceRatios;
	}


	public  LinkedHashMap<EnzymeName, Enzyme> getEnzymeHash(){
		return this.enzymesHash;
	}
	
	public  LinkedHashMap<MPathwayName, ArrayList<Enzyme>> getMetPathwaysHash(){
		return this.metPathwaysHash;
	}
	

}
