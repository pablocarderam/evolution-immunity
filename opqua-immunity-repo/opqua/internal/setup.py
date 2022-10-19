
"""Contains class Intervention."""

class Setup(object):
    """Class defines a setup with population parameters."""

    def __init__(
            self,
            id,
            num_loci, possible_alleles,
            fitnessHost, contactHost, receiveContactHost, mortalityHost,
            natalityHost, recoveryHost, migrationHost,
            populationContactHost, receivePopulationContactHost,
            mutationHost, recombinationHost,
            immunizationHost, deimmunizationHost,
            fitnessVector, contactVector, receiveContactVector, mortalityVector,
            natalityVector,recoveryVector, migrationVector,
            populationContactVector, receivePopulationContactVector,
            mutationVector, recombinationVector,
            immunizationVector, deimmunizationVector,
            contact_rate_host_vector,
            transmission_efficiency_host_vector,
            transmission_efficiency_vector_host,
            contact_rate_host_host,
            transmission_efficiency_host_host,
            mean_inoculum_host, mean_inoculum_vector,
            recovery_rate_host, recovery_rate_vector,
            mortality_rate_host,mortality_rate_vector,
            recombine_in_host, recombine_in_vector,
            num_crossover_host, num_crossover_vector,
            mutate_in_host, mutate_in_vector, death_rate_host,death_rate_vector,
            birth_rate_host, birth_rate_vector,
            vertical_transmission_host, vertical_transmission_vector,
            inherit_immunity_host, inherit_immunity_vector,
            immunity_acquisition_rate_host, immunity_acquisition_rate_vector,
            immunity_loss_rate_host, immunity_loss_rate_vector,
            immunityWeightsHost, immunityWeightsVector):
        """Create a new Setup.

        Arguments:
        id -- key of the Setup inside model dictionary (String)
        num_loci -- length of each pathogen genome string (int > 0)
        possible_alleles -- set of possible characters in all genome string, or
            at each position in genome string (String or list of Strings with
            num_loci elements)
        fitnessHost -- function that evaluates relative fitness in head-to-head
            competition for different genomes within the same host
            (function object, takes a String argument and returns a number >= 0)
        contactHost -- function that returns coefficient modifying probability
            of a given host being chosen to be the infector in a contact event,
            based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        receiveContactHost -- function that returns coefficient modifying
            probability of a given host being chosen to be the infected in
            a contact event, based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        mortalityHost -- function that returns coefficient modifying death rate
            for a given host, based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        natalityHost -- function that returns coefficient modifying birth rate
            for a given host, based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        recoveryHost -- function that returns coefficient modifying recovery
            rate for a given host based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        migrationHost -- function that returns coefficient modifying migration
            rate for a given host based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        populationContactHost -- function that returns coefficient modifying
            population contact rate for a given host based on genome sequence of
            pathogen
            (function object, takes a String argument and returns a number 0-1)
        mutationHost -- function that returns coefficient modifying mutation
            rate for a given host based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        recombinationHost -- function that returns coefficient modifying
            recombination rate for a given host based on genome sequence of
            pathogen
            (function object, takes a String argument and returns a number 0-1)
        immunizationHost -- function that returns coefficient modifying
            immunization rate for a given host based on genome sequence of
            pathogen
            (function object, takes a String argument and returns a number 0-1)
        deimmunizationHost -- function that returns coefficient modifying
            deimmunization rate for a given host based on genome sequence of
            pathogen
            (function object, takes a String argument and returns a number 0-1)
        fitnessVector -- function that evaluates relative fitness in head-to-
            head competition for different genomes within the same vector
            (function object, takes a String argument and returns a number >= 0)
        contactVector -- function that returns coefficient modifying probability
            of a given vector being chosen to be the infector in a contact
            event, based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        receiveContactVector -- function that returns coefficient modifying
            probability of a given vector being chosen to be the infected in
            a contact event, based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        mortalityVector -- function that returns coefficient modifying death
            rate for a given vector, based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        natalityVector -- function that returns coefficient modifying birth rate
            for a given vector, based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        recoveryVector -- function that returns coefficient modifying recovery
            rate for a given vector based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        migrationVector -- function that returns coefficient modifying migration
            rate for a given vector based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        populationContactVector -- function that returns coefficient modifying
            population contact rate for a given vector based on genome sequence
            of pathogen
            (function object, takes a String argument and returns a number 0-1)
        mutationVector -- function that returns coefficient modifying mutation
            rate for a given vector based on genome sequence of pathogen
            (function object, takes a String argument and returns a number 0-1)
        recombinationVector -- function that returns coefficient modifying
            recombination rate for a given vector based on genome sequence of
            pathogen
            (function object, takes a String argument and returns a number 0-1)
        immunizationVector -- function that returns coefficient modifying
            immunization rate for a given vector based on genome sequence of
            pathogen
            (function object, takes a String argument and returns a number 0-1)
        deimmunizationVector -- function that returns coefficient modifying
            deimmunization rate for a given vector based on genome sequence of
            pathogen
            (function object, takes a String argument and returns a number 0-1)
        contact_rate_host_vector -- rate of host-vector contact events, not
            necessarily transmission, assumes constant population density;
            evts/time (number >= 0)
        transmission_efficiency_host_vector -- fraction of host-vector contacts
            that result in successful transmission
        transmission_efficiency_vector_host -- fraction of vector-host contacts
            that result in successful transmission
        contact_rate_host_host -- rate of host-host contact events, not
            necessarily transmission, assumes constant population density;
            evts/time (number >= 0)
        transmission_efficiency_host_host -- fraction of host-host contacts
                that result in successful transmission
        mean_inoculum_host -- mean number of pathogens that are transmitted from
            a vector or host into a new host during a contact event (int >= 0)
        mean_inoculum_vector -- mean number of pathogens that are transmitted
            from a host to a vector during a contact event (int >= 0)
        recovery_rate_host -- rate at which hosts clear all pathogens;
            1/time (number >= 0)
        recovery_rate_vector -- rate at which vectors clear all pathogens
            1/time (number >= 0)
        recovery_rate_vector -- rate at which vectors clear all pathogens
            1/time (number >= 0)
        mortality_rate_host -- rate at which infected hosts die from disease
            (number 0-1)
        mortality_rate_vector -- rate at which infected vectors die from
            disease (number 0-1)
        recombine_in_host -- rate at which recombination occurs in host;
            evts/time (number >= 0)
        recombine_in_vector -- rate at which recombination occurs in vector;
            evts/time (number >= 0)
        num_crossover_host -- mean of a Poisson distribution modeling the number
            of crossover events of host recombination events (number >= 0)
        num_crossover_vector -- mean of a Poisson distribution modeling the
            number of crossover events of vector recombination events
            (number >= 0)
        mutate_in_host -- rate at which mutation occurs in host; evts/time
            (number >= 0)
        mutate_in_vector -- rate at which mutation occurs in vector; evts/time
            (number >= 0)
        death_rate_host -- natural host death rate; 1/time (number >= 0)
        death_rate_vector -- natural vector death rate; 1/time (number >= 0)
        birth_rate_host -- infected host birth rate; 1/time (number >= 0)
        birth_rate_vector -- infected vector birth rate; 1/time (number >= 0)
        vertical_transmission_host -- probability that a host is infected by its
            parent at birth (number 0-1)
        vertical_transmission_vector -- probability that a vector is infected by
            its parent at birth (number 0-1)
        inherit_immunity_host -- probability that a host inherits all
            immunity sequences from its parent (number 0-1)
        inherit_immunity_vector -- probability that a vector inherits all
            immunity sequences from its parent (number 0-1)
        immunity_acquisition_rate_host -- rate at which infected hosts acquire
            immunity to a pathogen infecting them; 1/time (number >= 0)
        immunity_acquisition_rate_vector -- rate at which infected vectors
            acquire immunity to a pathogen infecting them; 1/time (number >= 0)
        immunity_loss_rate_host -- rate at which infected hosts lose
            immunity to a pathogen infecting them; 1/time (number >= 0)
        immunity_loss_rate_vector -- rate at which infected vectors
            lose immunity to a pathogen infecting them; 1/time (number >= 0)
        immunityWeightsHost -- function that returns coefficient modifying
            immunity for a given host based on genome sequence of pathogen and a
            given immunity sequence
            (function object, takes two String arguments and returns a number)
        immunityWeightsVector -- function that returns coefficient modifying
            immunity for a given vector based on genome sequence of pathogen and
            a given immunity sequence
            (function object, takes two String arguments and returns a number)
        """

        super(Setup, self).__init__()

        self.id = id

        self.num_loci = num_loci
        if isinstance(possible_alleles, list):
            self.possible_alleles = possible_alleles
        else:
            self.possible_alleles = [possible_alleles] * self.num_loci
                # possible_alleles must be a list with all available alleles for
                # each position

        self.fitnessHost = fitnessHost
        self.contactHost = contactHost
        self.receiveContactHost = receiveContactHost
        self.mortalityHost = mortalityHost
        self.natalityHost = natalityHost
        self.recoveryHost = recoveryHost
        self.migrationHost = migrationHost
        self.populationContactHost = populationContactHost
        self.receivePopulationContactHost = receivePopulationContactHost
        self.mutationHost = mutationHost
        self.recombinationHost = recombinationHost
        self.immunizationHost = immunizationHost
        self.deimmunizationHost = deimmunizationHost

        self.fitnessVector = fitnessVector
        self.contactVector = contactVector
        self.receiveContactVector = receiveContactVector
        self.mortalityVector = mortalityVector
        self.natalityVector = natalityVector
        self.recoveryVector = recoveryVector
        self.migrationVector = migrationVector
        self.populationContactVector = populationContactVector
        self.receivePopulationContactVector = receivePopulationContactVector
        self.mutationVector = mutationVector
        self.recombinationVector = recombinationVector
        self.immunizationVector = immunizationVector
        self.deimmunizationVector = deimmunizationVector

        self.contact_rate_host_vector = contact_rate_host_vector
        self.contact_rate_host_host = contact_rate_host_host
            # contact rates assumes scaling area--large populations are equally
            # dense as small ones, so contact is constant with both host and
            # vector populations. If you don't want this to happen, modify the
            # population's contact rate accordingly.
            # Examines contacts between infected hosts and all hosts
        self.transmission_efficiency_host_vector = transmission_efficiency_host_vector
        self.transmission_efficiency_vector_host = transmission_efficiency_vector_host
        self.transmission_efficiency_host_host = transmission_efficiency_host_host
        self.mean_inoculum_host = mean_inoculum_host
        self.mean_inoculum_vector = mean_inoculum_vector
        self.recovery_rate_host = recovery_rate_host
        self.recovery_rate_vector = recovery_rate_vector
        self.mortality_rate_host = mortality_rate_host
        self.mortality_rate_vector = mortality_rate_vector

        self.recombine_in_host = recombine_in_host
        self.recombine_in_vector = recombine_in_vector
        self.num_crossover_host = num_crossover_host
        self.num_crossover_vector = num_crossover_vector
        self.mutate_in_host = mutate_in_host
        self.mutate_in_vector = mutate_in_vector

        self.death_rate_host = death_rate_host
        self.death_rate_vector = death_rate_vector
        self.birth_rate_host = birth_rate_host
        self.birth_rate_vector = birth_rate_vector

        self.vertical_transmission_host = vertical_transmission_host
        self.vertical_transmission_vector = vertical_transmission_vector
        self.inherit_immunity_host = inherit_immunity_host
        self.inherit_immunity_vector = inherit_immunity_vector

        self.immunity_acquisition_rate_host = immunity_acquisition_rate_host
        self.immunity_acquisition_rate_vector = immunity_acquisition_rate_vector
        self.immunity_loss_rate_host = immunity_loss_rate_host
        self.immunity_loss_rate_vector = immunity_loss_rate_vector

        self.immunityWeightsHost = immunityWeightsHost
        self.immunityWeightsVector = immunityWeightsVector

        # TODO: add coefficient weights for each genome position
