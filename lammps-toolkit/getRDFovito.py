#!/usr/bin/env python3
import argparse

from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier,TimeAveragingModifier

parser = argparse.ArgumentParser(description="calculate RDF by Ovito python")
parser.add_argument("-i", "--filename", default="XDATCAR", help="Path to the XDATCAR file (default: 'XDATCAR'). If you use vasp-format XDATCAR, the name can only be XDATCAR. ")
parser.add_argument("-c", "--cutoff", default=5.0, type=float, help="cutoff radius")
parser.add_argument("-n", "--num_of_bins", default=100, type=int, help="number of bins")
parser.add_argument("-o", "--outputname", default='rdf.txt', type=str, help="output filename")
parser.add_argument("-p", "--partial", default=False, type=bool, help="whether to get partial element or not")
args = parser.parse_args()

pipeline = import_file(args.filename)
print("Number of MD frames:", pipeline.num_frames)

# Print the list of input particle types.
# They are represented by ParticleType objects attached to the 'Particle Type' particle property.
for t in pipeline.compute().particles.particle_types.types:
    print("Type %i: %s" % (t.id, t.name))

# Insert the RDF calculation modifier into the pipeline:
pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=args.cutoff, number_of_bins=args.num_of_bins, partial=args.partial))

# Insert the time-averaging modifier into the pipeline, which accumulates
# the instantaneous DataTable produced by the previous modifier and computes a mean histogram.
pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))

# Data export method 1: Convert to NumPy array and write data to a text file:
export_file(pipeline,args.outputname, "txt/table", key="coordination-rdf[average]")

#modifier = CoordinationAnalysisModifier(cutoff=args.cutoff, number_of_bins=args.num_of_bins, partial=True)
#pipeline.modifiers.append(modifier)
#pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
#export_file(pipeline,args.outputname,"txt/table",key="coordination-rdf[average]")
