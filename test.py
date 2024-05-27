from radical.radical import Radical
a=Radical("sk",export_uncal_pdf=True, export_resage_pdf=True,threshold=1e-14, mixture_pdf=True)
a.radcal()