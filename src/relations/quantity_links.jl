# relations/quantity_links.jl — the relation → quantity map.
#
# Declares, for each relation whose subject is a named quantity, which
# `AbstractQuantity` TYPES it directly constrains (parametric families use
# the UnionAll, so `Susceptibility` matches any `Susceptibility{I}`).  This
# is the machine link that turns the registry into a queryable web:
# `relations_constraining(q)` reads it in reverse.  Included after the
# physics submodules are re-exported, so every relation + quantity name
# resolves at the top level.  Relations that constrain parameters or
# exponents rather than named quantities (scaling laws, Maxwell relations,
# the uncertainty relation) keep the default `quantities(rel) == ()`.

# ── Statistical mechanics ──
quantities(::MagnetizationResponse) = (Magnetization,)
quantities(::SusceptibilityResponse) = (Susceptibility,)
quantities(::SusceptibilityFDT) = (Susceptibility, Magnetization)
quantities(::SusceptibilityPositivity) = (Susceptibility,)
quantities(::StructureFactorSusceptibility) = (Susceptibility, StaticStructureFactor)
quantities(::SpecificHeatFDT) = (SpecificHeat, Energy)
quantities(::SpecificHeatFromEntropy) = (SpecificHeat, ThermalEntropy)
quantities(::SpecificHeatPositivity) = (SpecificHeat,)
quantities(::HeatCapacityDifference) = (SpecificHeat, IsobaricSpecificHeat)
quantities(::CompressibilityPositivity) = (IsothermalCompressibility,)
quantities(::EntropyResponse) = (ThermalEntropy, FreeEnergy)
quantities(::GibbsHelmholtz) = (Energy, FreeEnergy)
quantities(::FreeEnergyFromZ) = (FreeEnergy, PartitionFunction)
quantities(::FreeEnergyLegendre) = (FreeEnergy, Energy, ThermalEntropy)
quantities(::ClausiusClapeyron) = (LatentHeat,)
quantities(::GibbsDuhem) = (ThermalEntropy, Volume, ParticleNumber)
quantities(::MicrocanonicalTemperature) = (ThermalEntropy, Energy)
quantities(::CanonicalTPQ) = (PartitionFunction,)

# ── Correlations / Green's functions ──
quantities(::Dyson) = (RetardedGreensFunction, SelfEnergy)
quantities(::SpectralFromGreens) = (SpectralFunction, RetardedGreensFunction)
quantities(::SpectralSumRule) = (SpectralFunction,)
quantities(::DetailedBalance) = (DynamicalStructureFactor,)
quantities(::DynamicalFDT) = (DynamicalStructureFactor, DynamicalSusceptibility)
quantities(::CorrelationLengthGap) = (CorrelationLength, MassGap)
quantities(::NMRExponent) = (NMRSpinRelaxationRate,)
quantities(::FiniteSizeGap) = (MassGap,)
function quantities(::StaticFromDynamicalStructureFactor)
    return (StaticStructureFactor, DynamicalStructureFactor)
end

# ── Transport ──
quantities(::WiedemannFranz) = (ThermalConductivity, Conductivity)
quantities(::RighiLeduc) = (ThermalConductivity,)
quantities(::MottFormula) = (Thermopower, Conductivity)
quantities(::KelvinRelation) = (PeltierCoefficient, Thermopower)
quantities(::OpticalSumRule) = (DrudeWeight,)
quantities(::CurrentNoiseFDT) = (CurrentNoise,)
quantities(::ThermoelectricFigureOfMerit) = (Thermopower, Conductivity, ThermalConductivity)
quantities(::LongitudinalResistivity) = (Resistivity, Conductivity)
quantities(::HallResistivity) = (Resistivity, Conductivity)
quantities(::VonKlitzing) = (Resistivity, FillingFactor)
quantities(::MobilityConductivity) = (Conductivity, Mobility, CarrierDensity)
quantities(::EinsteinRelation) = (Mobility, DiffusionConstant)
quantities(::HallAngle) = (Conductivity,)
quantities(::CyclotronFrequency) = (EffectiveMass, MagneticFluxDensity)
quantities(::PowerFactor) = (Thermopower, Conductivity)
quantities(::NernstCoefficient) = (Thermopower, MagneticFluxDensity)

# ── Quantum information & entanglement ──
quantities(::RenyiTwoPurity) = (RenyiEntropy, Purity)
quantities(::RenyiEntropyMoment) = (RenyiEntropy,)
quantities(::TsallisEntropyMoment) = (TsallisEntropy,)
quantities(::MutualInformationDefinition) = (MutualInformation,)
quantities(::ConditionalEntropyDefinition) = (ConditionalEntropy,)
quantities(::MeasurementEntropyIncrease) = (MeasurementEntropy,)
quantities(::MeasurementEntropyRelative) = (MeasurementEntropy, RelativeEntropy)
quantities(::MarkovEntropyDefinition) = (MarkovEntropy,)
quantities(::ConcurrenceTangle) = (Concurrence, Tangle)
quantities(::Monogamy) = (Tangle,)
quantities(::ThreeTangleDefinition) = (ThreeTangle, Tangle)
quantities(::TripartiteInformationDefinition) = (TripartiteInformation,)
quantities(::KitaevPreskillTEE) = (TopologicalEntanglementEntropy,)
quantities(::EntropyNonNegativity) = (VonNeumannEntropy,)
quantities(::MaxEntropyBound) = (VonNeumannEntropy,)
quantities(::Subadditivity) = (VonNeumannEntropy,)
quantities(::ArakiLieb) = (VonNeumannEntropy,)
quantities(::StrongSubadditivity) = (VonNeumannEntropy,)
quantities(::RenyiMonotonicity) = (RenyiEntropy,)
quantities(::RelativeEntropyNonNegativity) = (RelativeEntropy,)
quantities(::CFTEntanglementSlope) = (VonNeumannEntropy,)

# ── Quantum-mechanical foundations ──
quantities(::VirialTheorem) = (KineticEnergy, PotentialEnergy)
quantities(::EnergyVarianceEigenstate) = (EnergyVariance,)

# ── Topology ──
quantities(::TKNN) = (Conductivity,)
