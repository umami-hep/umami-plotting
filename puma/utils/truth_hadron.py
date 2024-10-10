import numpy as np

def MaskTracks(my_data, n_jets, n_tracks):

    n_real_tracks = np.repeat(my_data["jets"]["n_tracks"], n_tracks).reshape(n_jets, n_tracks) # This is needed because jets have a different format than tracks
    track_indices = np.tile(
         np.arange(0,n_tracks,dtype=np.int32),
         n_jets,
    ).reshape(n_jets, n_tracks)
    
    track_mask = np.where(track_indices < n_real_tracks, 1, 0)
    
    return track_mask, n_real_tracks


def ProcessTruthHadrons(my_data, n_jets):

    hadron = my_data["truth_hadrons"][:] 
    max_n_hadrons = hadron["ftagTruthParentBarcode"].shape[1]
    
    print("Max number of truth hadrons in your sample: ", max_n_hadrons)

    # Make your own loop
    child = []
    parent = []

    for i in range(0, n_jets):
        child.append(np.isin(hadron["ftagTruthParentBarcode"][i], hadron["barcode"][i]))
        parent.append(np.isin(hadron["barcode"][i], hadron["ftagTruthParentBarcode"][i]))
            
    child = np.where(hadron["ftagTruthParentBarcode"] == -1, False,np.array(child))
    parent = np.where(hadron["barcode"] == -1, False, np.array(parent))
    
    real_hadron = np.where(hadron["barcode"] == -1, False, True)
    one_hadron = np.where(np.sum(real_hadron, axis =1) ==1, True, False)

    label_hadron_index = np.where(np.sum(parent) > 0, np.argmax(parent, axis=1), np.nan)
    label_hadron_index = np.where(one_hadron == True, 0, label_hadron_index)

    child_hadron_index = np.where(np.sum(child) > 0 , np.argmax(child, axis=1), np.nan) # argmax takes the first argument from the list
    child_hadron_index = np.where(one_hadron == True, np.nan, child_hadron_index) # make sure that if there is only one hadron, the argmax does not select 0
    
    child2_hadron_index = np.full(child_hadron_index.shape, np.nan)

    for i in range (0, n_jets):
        if (np.sum(child[i]) >= 2): 
            for x in (np.argwhere(child[i] == 1)):
                if x != child_hadron_index[i]:
                    child2_hadron_index[i] = x
                else:
                    continue
            
    hadron_indices = np.stack( (label_hadron_index, child_hadron_index,  child2_hadron_index), axis = 1)

    print("Total number of jets = ", n_jets)
    n_one = np.sum(one_hadron)
    
    family = np.where(hadron["barcode"] == -1, False, child | parent)
    unrelated = ((np.sum(real_hadron, axis =1) > 1) & (np.sum(family, axis =1) <= 1))

    good_jets = ((one_hadron) | (np.sum(family, axis =1) > 1))

    n_decay_chain = np.sum(np.where(np.sum(family, axis =1) > 1, True, False))
    n_unrelated = np.sum(unrelated)
    print(" -- Number of jets multiple unrelated hadrons = ",n_unrelated)

    print(" -- Number of jets with only 1 SV = ",n_one)

    print(" -- Number of jets with a decay chain = ",n_decay_chain)
    if ( n_one + n_decay_chain + n_unrelated != n_jets): 
        missing_jets = hadron[np.invert(good_jets) & np.invert(unrelated)]        
        print("The sum of the 3 types of jets does not add up to the total number of jets!")
        if (np.sum(np.where(missing_jets["pdgId"] != -1, 1, 0)) == 0): 
            if (np.sum(np.where(missing_jets["flavour"] != -1, 1, 0)) == 0):
                print("This is because for ", len(missing_jets), " non-HF jets the pdgID of all truth hadrons is -1")

        else:
            print("The reason is not understood for ", np.sum(np.where(missing_jets["pdgId"] != -1, 1, 0)), " jets and you should DEBUG this")

    print("   --> # jets with 2 SVs: = ",np.sum(np.where(np.sum(family, axis =1) == 2, True, False)))
    print("   --> # jets with 3+ SVs: = ",np.sum(np.where(np.sum(family, axis =1) > 2, True, False)))
    print("   --> # jets with 4+ SVs: = ",np.sum(np.where(np.sum(family, axis =1) > 3, True, False)))

    print(np.sum(good_jets), " good jets (= the jets with more than one unrelated hadron and light jets with no truth hadrons have been dropped )")

    
    return good_jets, hadron_indices, parent, child, one_hadron, unrelated


def AssociateTracksToHadron(my_data, good_jets, drop_bad_jets = True, debug=False):
    
    n_jets, n_tracks = my_data["tracks"].shape
    #track = my_data["tracks"][:, 0:n_tracks]
    
    # Start by getting a mask of the real tracks
    # Get real tracks
    track_mask, n_real_tracks = MaskTracks(my_data, n_jets, n_tracks)
    n_hadrons = my_data["truth_hadrons"].shape[1]
    
    if drop_bad_jets:
        my_good_jets = np.repeat(good_jets, n_tracks).reshape(n_jets, n_tracks) # This is needed because jets have a different format than tracks
        jet_track_mask = np.where(my_good_jets == True, 1, 0)
        track_mask = track_mask & jet_track_mask

    inclusive_vertex = []
    exclusive_vertex = []
    n_tracks_inclusive_vertex = []
    n_tracks_exclusive_vertex = []
    good_hadron_track_association = []

    hadron_index = []

    dummy = np.zeros(n_tracks).astype(int)

    for i in range(0, n_jets):
        if drop_bad_jets:
            if good_jets[i] == False: continue


        if debug: print("Track Parent Barcodes ", my_data["tracks"]["ftagTruthParentBarcode"][i] )
        positive_track_barcodes = np.where(my_data["tracks"]["ftagTruthParentBarcode"][i] < 0, np.nan, my_data["tracks"]["ftagTruthParentBarcode"][i])
        if debug: print("Hadron Barcodes ", my_data["truth_hadrons"]["barcode"][i])
        inclusive_tracks_to_hadron = np.where(track_mask[i] == 0, 0, np.isin(positive_track_barcodes,  my_data["truth_hadrons"]["barcode"][i])).astype(int)        
        n_tracks_inclusive_vertex.append(np.sum(inclusive_tracks_to_hadron))
        if np.sum(inclusive_tracks_to_hadron) <= 1: 
            inclusive_vertex.append(dummy)
        else:
            inclusive_vertex.append(inclusive_tracks_to_hadron)        
        
        tmp_exclusive_list = []
        tmp_n_tracks = []
        for j in range(0, n_hadrons):
            tmp_exclusive = np.where(track_mask[i] == 0, 0, np.isin(positive_track_barcodes,  my_data["truth_hadrons"]["barcode"][i][j])).astype(int)        
            if np.sum(tmp_exclusive) <= 1:    
                
                tmp_exclusive_list.append(dummy)
            else:
                tmp_exclusive_list.append(tmp_exclusive)
            tmp_n_tracks.append(int(np.sum(tmp_exclusive_list[-1])))
        n_tracks_exclusive_vertex.append(tmp_n_tracks)
    
        if np.sum(n_tracks_exclusive_vertex[-1]) > 0:        
            hadron_index.append(np.argmax(n_tracks_exclusive_vertex[-1]))
            good_hadron_track_association.append(1)
        else:
            hadron_index.append(-99) # None of the hadron indices have 2 associated trakcks
            good_hadron_track_association.append(0)


        exclusive_vertex.append(np.array(tmp_exclusive_list)) 

        if debug: print("INCLUSIVE Track to Hadron Association ", inclusive_vertex[-1])
        if debug: print("EXCLUSIVE Track to Hadron Association ", exclusive_vertex[-1])
        if debug: print("Index for Hadron with highest amoutn of associated tracks ", hadron_index[-1])
        if debug: print("Number of tracks associated to the each hadron ", n_tracks_per_hadron[-1])
        
    inclsuive_vertex = np.array(inclusive_vertex)
    exclusive_vertex = np.array(exclusive_vertex)
    hadron_index = np.array(hadron_index)

    
    return inclusive_vertex, exclusive_vertex, hadron_index, n_tracks_inclusive_vertex, n_tracks_exclusive_vertex, track_mask, jet_track_mask, good_hadron_track_association

def SelectHadronMostTracks(truth_hadrons, hadron_index):
    invalid_jet_mask = np.where(hadron_index == -99, 1, 0).astype(int)
    
    # remove -99 so that you can get the index, but you will need to add the mask later
    tmp_hadron_index = np.where(invalid_jet_mask == 1, 0, hadron_index)  

    # select hadron with most tracks
    truth_hadron_most_tracks = truth_hadrons[np.arange(truth_hadrons.shape[0]),np.array(tmp_hadron_index).astype(int)]

    # mask invalid jets
    truth_hadron_most_tracks = truth_hadron_most_tracks[~invalid_jet_mask]

    return truth_hadron_most_tracks
